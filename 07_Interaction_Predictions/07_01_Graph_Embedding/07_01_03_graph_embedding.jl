using Pkg
Pkg.activate(".")
using SpeciesInteractionNetworks
using CSV
using DataFrames
using StatsBase
using Pipe
using Distributions
using PhyloNetworks #Pkg.add(name="PhyloNetworks", version="0.16.4")
using Random
using JLD2



#Read in the data
interaction_df = @pipe CSV.read("Data_Clean/Interactions/ints_final_clean_slim.csv", DataFrame) |>
                    filter(row -> row.interactionTypeName .!= "hasHost", _) |> #filter out crops for now
                    filter(row -> length(split(row.sourceTaxonName_harm)) >= 2, _) # filter out pollinators that is at the genus level only

#check the interactionTypeName (no host!)
unique(interaction_df.interactionTypeName)

#we need to extract the first and second word from sourceTaxonName_harm
interaction_df.genus_species=[join(split(name)[1:2], " ") for name in interaction_df.sourceTaxonName_harm]

checklist_df = @pipe CSV.read("Data_Clean/Species_Checklists/checklist_cleaned.csv", DataFrame) 


#Filter the checklist to only include the targeted group
group_name = "moths" #change this to the targeted group
group_checklist = checklist_df[checklist_df.taxon .∈ Ref([group_name]), :]

#Filter the interaction data to only include the targeted group based the genera present in the checklist
group_genera = unique(group_checklist.genus)
group_interaction_df = @pipe interaction_df |>
    filter(row -> row.sourceTaxonGenusName ∈ group_genera, _)

#List of the unique species in the interaction data for the targeted group
group_species = unique(group_interaction_df.genus_species)

#Now we need the plant genera in the interaction data
plant_genera = unique(group_interaction_df.targetTaxonGenusName)
#Filter the checklist to only include the plant genera present in the interaction data
plant_checklist = @pipe checklist_df |>
    filter(row -> row.genus ∈ plant_genera, _)

#Read in the phylogenetic tree for the targeted group
pollinator_tree = readTopology("Data_Clean/Phylogeny/$(group_name)_tree.nwk")

# A function to perform ancestral state reconstruction
function ace_distribution(τ, ϕ)
    astate = ancestralStateReconstruction(τ, ϕ)
    ustate = predint(astate, level=0.3)
    pstate = predint(astate, level=0.95)
    Du = Uniform.(ustate[:,1].-10eps(), ustate[:,2].+10eps())
    σ = (pstate[:,2] .- pstate[:,1])./3.92
    μ = mean(pstate,dims = 2)[:,1]
    Dn = Normal.(μ, σ)
    return (Du, Dn)
end


#We will compare between 1.0 and 0.8 variance explained
#the 10-fold CV for getting the average threshold
n_folds = 10
var_explained = [0.8, 1.0] #if it is 1.0 then there will be no in-sample learning, only out-of-sample prediction, but if it is less than 1.0 then there will be some in-sample learning


fold_metrics_df = DataFrame()

for v in var_explained
    println("Variance explained: ", v)

    for fold in 1:n_folds
        println("Fold $fold/$n_folds")
        
        # Hold out 10% of the pollinator species
        Random.seed!(fold)
        n_species = length(group_species)
        n_holdout = max(2, Int(floor(0.1 * n_species))) #ensure at least two species are held out
        
        shuffled_species = shuffle(group_species)
        holdout_species = shuffled_species[1:n_holdout]
        train_species = shuffled_species[n_holdout+1:end]
        
        println("Held-out genera: ", join(holdout_species, ", "))
        
        # Filter interactions to training species only
        group_interaction_train = @pipe group_interaction_df |>
            filter(row -> !(row.genus_species ∈ holdout_species), _)
        
        # Build training network
        pollinator_in_interaction_train = String.(unique(group_interaction_train.genus_species))

        plant_in_interaction_train = String.(unique(group_interaction_train.targetTaxonGenusName))
        
        nodes = Bipartite(pollinator_in_interaction_train, plant_in_interaction_train)
        
        #look at which plant genera are missing from the training data
        missing_plant_genera = setdiff(plant_genera, unique(plant_in_interaction_train))

        #Set up the training interaction matrix
        interaction_matrix_train = zeros(Bool, length(pollinator_in_interaction_train), length(plant_in_interaction_train))
        
        for row in eachrow(group_interaction_train)
            pollinator_idx = findfirst(==(row.genus_species), pollinator_in_interaction_train)
            plant_idx = findfirst(==(row.targetTaxonGenusName), plant_in_interaction_train)
            
            if pollinator_idx !== nothing && plant_idx !== nothing
                interaction_matrix_train[pollinator_idx, plant_idx] = true
            end
        end
        
        edges_train = Binary(interaction_matrix_train)
        M_train = SpeciesInteractionNetwork(nodes, edges_train)
        
        # Determine k from variance threshold
        rnk = rank(M_train)
        eig = svd(M_train).S[1:rnk]

        neig = eig ./ sum(eig)
        cumvar = cumsum(neig)
        k = findfirst(round.(cumvar, digits=10) .>= v)

        # Embed training network
        L, R = tsvd(M_train, k)
        # Now use phylogeny to predict latent variables for all the held-out species
        pollinator_leaves = tipLabels(pollinator_tree)
        pollinator_pool = DataFrame(; tipNames=unique(replace.(holdout_species, " " => "_"))) #the held-out species 
        pollinator_traitframe = DataFrame(; tipNames=pollinator_leaves)
        
        leftnames = "L" .* string.(1:size(L, 2))
        traits_L = DataFrame(L, Symbol.(leftnames))
        traits_L[!, "tipNames"] = replace.(pollinator_in_interaction_train, " " => "_")
        
        pollinator_traits = leftjoin(pollinator_traitframe, traits_L; on=:tipNames)
        
        # Reconstruct latent variables for all species (including held-out)
        pollinator_imputedtraits = DataFrame(;
            tipNames=[pollinator_leaves; fill(missing, pollinator_tree.numNodes - pollinator_tree.numTaxa)]

        )
        
        # Phylogenetic reconstruction for each dimension
        for coord in 1:k
            trait_col = pollinator_traits[!, ["L$(coord)", "tipNames"]]
            Du, Dn = ace_distribution(trait_col, pollinator_tree)
            pollinator_imputedtraits[!, "L$(coord)_Normal"] = Dn
        end

        pollinator_rec = innerjoin(dropmissing(pollinator_imputedtraits), pollinator_pool; on=:tipNames)#some held-out species may not be in the phylogeny

        #Print out the held-out species that are not in the phylogeny
        missing_in_phylo = setdiff(holdout_species, replace.(pollinator_rec.tipNames, "_" => " "))
        if !isempty(missing_in_phylo)
            println("The following held-out genera are not in the phylogeny and will be excluded from evaluation: ", join(missing_in_phylo, ", "))
        end
                    
        # Sample from phylogeny-inferred distributions
        draws = 20000
        ℒn = Array(pollinator_rec[!, leftnames .* "_Normal"])
        
        # Predict interactions for held-out species
        # Average the latent variable over the 2000 draws
        L_sum = zeros(size(ℒn))
        for i in 1:draws
            L_sum .+= rand.(ℒn)
        end
        L_mean = L_sum ./ draws
        
        #calculate the dot product between L_mean and R
        predicted_dp = L_mean * R

        #Get the observed interaction natrix for the held-out species only
        interaction_matrix_validation = zeros(Bool, length(pollinator_rec.tipNames), length(plant_in_interaction_train))
        
        for row in eachrow(group_interaction_df)
            pollinator_idx = findfirst(==(row.genus_species), replace.(pollinator_rec.tipNames, "_" => " "))
            plant_idx = findfirst(==(row.targetTaxonGenusName), plant_in_interaction_train)
            
            if pollinator_idx !== nothing && plant_idx !== nothing
                interaction_matrix_validation[pollinator_idx, plant_idx] = true
            end
        end

        #Now we can calculate the performance metrics based on different threshold values
        thresholds = LinRange(extrema(predicted_dp)..., 500)
                    
        tp = zeros(Float64, length(thresholds))
        tn = similar(tp)
        fp = similar(tp)
        fn = similar(tp)
        
        A = interaction_matrix_validation
        # For each threshold, calculate performance metrics
        for (i, t) in enumerate(thresholds)
            PN = predicted_dp .>= t # Predicted network at threshold t
            tp[i] = sum((A) .& (PN)) / sum(A) # True positive rate
            tn[i] = sum((.!(A)) .& (.!(PN))) / sum(.!(A)) # True negative rate
            fp[i] = sum((.!(A)) .& (PN)) / sum(.!(A)) # False positive rate
            fn[i] = sum((A) .& (.!(PN))) / sum(A) # False negative rate
        end
        
        # Find optimal threshold using Youden's J
        Y = tp .+ tn .- 1
        maxY, posY = findmax(Y)
        threshold = thresholds[posY]
        println("Optimal threshold: ", threshold)
                                

        # Evaluate performance on held-out genera
        PN_final = predicted_dp .>= threshold
        tp_final = sum((A) .& (PN_final))
        tn_final = sum((.!(A)) .& (.!(PN_final)))
        fp_final = sum((.!(A)) .& (PN_final))
        fn_final = sum((A) .& (.!(PN_final)))



        #return a dataframe with the metrics
        fold_metric = DataFrame(
            fold = fold,
            variance_explained = v,
            threshold = threshold,
            recall = tp_final/(tp_final+fn_final),
            precision = tp_final/(tp_final+fp_final))    

        append!(fold_metrics_df, fold_metric) 
    
    end
    
end

#save the results to a csv file
CSV.write(
    "Data_Clean/Graph_Embedding/$(group_name)_cv_for_selecting_threshold.csv",
    fold_metrics_df
)


#Then now we can use the average threshold from the above to evaluate the performance on the interaction data, we will do a 3-fold cross validation here
avg_threshold_0_8 = mean(fold_metrics_df.threshold[fold_metrics_df.variance_explained .== 0.8])
avg_threshold_1_0 = mean(fold_metrics_df.threshold[fold_metrics_df.variance_explained .== 1.0])

n_folds = 3
performance_metrics_df = DataFrame()

for v in var_explained
    println("Variance explained: ", var)

    if v == 0.8
        avg_threshold = avg_threshold_0_8
    else
        avg_threshold = avg_threshold_1_0
    end

    for fold in 1:n_folds

        

        println("Fold $fold/$n_folds")
        
        # Hold out 10% of the pollinator species
        #setup a random seed_value
        seed_value = rand(1:10000)
        Random.seed!(seed_value)
        n_species = length(group_species)
        n_holdout = max(2, Int(floor(0.1 * n_species))) #ensure at least two species are held out
        
        shuffled_species = shuffle(group_species)
        holdout_species = shuffled_species[1:n_holdout]
        train_species = shuffled_species[n_holdout+1:end]
        
        println("Held-out genera: ", join(holdout_species, ", "))
        
        # Filter interactions to training species only
        group_interaction_train = @pipe group_interaction_df |>
            filter(row -> !(row.genus_species ∈ holdout_species), _)
        
        # Build training network
        pollinator_in_interaction_train = String.(unique(group_interaction_train.genus_species))

        plant_in_interaction_train = String.(unique(group_interaction_train.targetTaxonGenusName))
        
        nodes = Bipartite(pollinator_in_interaction_train, plant_in_interaction_train)
        
        #look at which plant genera are missing from the training data
        missing_plant_genera = setdiff(plant_genera, unique(plant_in_interaction_train))

        #Set up the training interaction matrix
        interaction_matrix_train = zeros(Bool, length(pollinator_in_interaction_train), length(plant_in_interaction_train))
        
        for row in eachrow(group_interaction_train)
            pollinator_idx = findfirst(==(row.genus_species), pollinator_in_interaction_train)
            plant_idx = findfirst(==(row.targetTaxonGenusName), plant_in_interaction_train)
            
            if pollinator_idx !== nothing && plant_idx !== nothing
                interaction_matrix_train[pollinator_idx, plant_idx] = true
            end
        end
        
        edges_train = Binary(interaction_matrix_train)
        M_train = SpeciesInteractionNetwork(nodes, edges_train)
        
        # Determine k from variance threshold
        rnk = rank(M_train)
        eig = svd(M_train).S[1:rnk]

        neig = eig ./ sum(eig)
        cumvar = cumsum(neig)
        k = findfirst(round.(cumvar, digits=10) .>= v)

        # Embed training network
        L, R = tsvd(M_train, k)
        # Now use phylogeny to predict latent variables for all the held-out species
        pollinator_leaves = tipLabels(pollinator_tree)
        pollinator_pool = DataFrame(; tipNames=unique(replace.(holdout_species, " " => "_"))) #the held-out species 
        pollinator_traitframe = DataFrame(; tipNames=pollinator_leaves)
        
        leftnames = "L" .* string.(1:size(L, 2))
        traits_L = DataFrame(L, Symbol.(leftnames))
        traits_L[!, "tipNames"] = replace.(pollinator_in_interaction_train, " " => "_")
        
        pollinator_traits = leftjoin(pollinator_traitframe, traits_L; on=:tipNames)
        
        # Reconstruct latent variables for all species (including held-out)
        pollinator_imputedtraits = DataFrame(;
            tipNames=[pollinator_leaves; fill(missing, pollinator_tree.numNodes - pollinator_tree.numTaxa)]
        )
        
        # Phylogenetic reconstruction for each dimension
        for coord in 1:k
            trait_col = pollinator_traits[!, ["L$(coord)", "tipNames"]]
            Du, Dn = ace_distribution(trait_col, pollinator_tree)
            pollinator_imputedtraits[!, "L$(coord)_Normal"] = Dn
        end

        pollinator_rec = innerjoin(dropmissing(pollinator_imputedtraits), pollinator_pool; on=:tipNames)#some held-out species may not be in the phylogeny

        #Print out the held-out species that are not in the phylogeny
        missing_in_phylo = setdiff(holdout_species, replace.(pollinator_rec.tipNames, "_" => " "))
        if !isempty(missing_in_phylo)
            println("The following held-out genera are not in the phylogeny and will be excluded from evaluation: ", join(missing_in_phylo, ", "))
        end
                    
        # Sample from phylogeny-inferred distributions
        draws = 20000
        ℒn = Array(pollinator_rec[!, leftnames .* "_Normal"])
        
        # Predict interactions for held-out species
        # Average the latent variable over the 2000 draws
        L_sum = zeros(size(ℒn))
        for i in 1:draws
            L_sum .+= rand.(ℒn)
        end
        L_mean = L_sum ./ draws

        
        #calculate the dot product between L_mean and R
        predicted_dp = L_mean * R

        #Get the observed interaction matrix for the held-out species only
        interaction_matrix_validation = zeros(Bool, length(pollinator_rec.tipNames), length(plant_in_interaction_train))
        
        for row in eachrow(group_interaction_df)
            pollinator_idx = findfirst(==(row.genus_species), replace.(pollinator_rec.tipNames, "_" => " "))
            plant_idx = findfirst(==(row.targetTaxonGenusName), plant_in_interaction_train)
            
            if pollinator_idx !== nothing && plant_idx !== nothing
                interaction_matrix_validation[pollinator_idx, plant_idx] = true
            end
        end

        A = interaction_matrix_validation

        # Evaluate performance on held-out genera using the average threshold
        PN_final = predicted_dp .>= avg_threshold
        tp_final = sum((A) .& (PN_final))
        tn_final = sum((.!(A)) .& (.!(PN_final)))
        fp_final = sum((.!(A)) .& (PN_final))
        fn_final = sum((A) .& (.!(PN_final)))

        #return a dataframe with the metrics
        performance_metric = DataFrame(
            fold = fold,
            variance_explained = v,
            threshold = avg_threshold,
            recall = tp_final/(tp_final+fn_final),
            precision = tp_final/(tp_final+fp_final))
        append!(performance_metrics_df, performance_metric)
    end
    
end
    
#save the results to a csv file
CSV.write(
    "Data_Clean/Graph_Embedding/$(group_name)_final_validation_performance.csv",
    performance_metrics_df
)



#Finally, we can now predict interactions for the species on the checklist using the two variance explained levels and the average threshold from the cross-validation above

for v in var_explained
    println("Variance explained: ", var)

    if v == 0.8
        avg_threshold = avg_threshold_0_8
    else
        avg_threshold = avg_threshold_1_0
    end

        
    # Build the full network
    pollinator_in_interaction = String.(unique(group_interaction_df.genus_species))

    plant_in_interaction = String.(unique(group_interaction_df.targetTaxonGenusName))
    
    nodes = Bipartite(pollinator_in_interaction, plant_in_interaction)
    
    #Set up the training interaction matrix
    interaction_matrix = zeros(Bool, length(pollinator_in_interaction), length(plant_in_interaction))
    
    for row in eachrow(group_interaction_df)
        pollinator_idx = findfirst(==(row.genus_species), pollinator_in_interaction)
        plant_idx = findfirst(==(row.targetTaxonGenusName), plant_in_interaction)
        
        if pollinator_idx !== nothing && plant_idx !== nothing
            interaction_matrix[pollinator_idx, plant_idx] = true
        end
    end
    
    edges = Binary(interaction_matrix)
    M = SpeciesInteractionNetwork(nodes, edges)
    
    # Determine k from variance threshold
    rnk = rank(M)
    eig = svd(M).S[1:rnk]

    neig = eig ./ sum(eig)
    cumvar = cumsum(neig)
    k = findfirst(round.(cumvar, digits=10) .>= v)

    # Embed training network
    L, R = tsvd(M, k)
    # Now use phylogeny to predict latent variables for all the held-out species
    pollinator_leaves = tipLabels(pollinator_tree)
    pollinator_pool = DataFrame(; tipNames=unique(replace.(group_checklist.genus_species, " " => "_"))) #the species on the checklist
    pollinator_traitframe = DataFrame(; tipNames=pollinator_leaves)
    
    leftnames = "L" .* string.(1:size(L, 2))
    traits_L = DataFrame(L, Symbol.(leftnames))
    traits_L[!, "tipNames"] = replace.(pollinator_in_interaction, " " => "_")
    
    pollinator_traits = leftjoin(pollinator_traitframe, traits_L; on=:tipNames)
    
    # Reconstruct latent variables for all species (including held-out)
    pollinator_imputedtraits = DataFrame(;
        tipNames=[pollinator_leaves; fill(missing, pollinator_tree.numNodes - pollinator_tree.numTaxa)]
    )
    
    # Phylogenetic reconstruction for each dimension
    for coord in 1:k
        trait_col = pollinator_traits[!, ["L$(coord)", "tipNames"]]
        Du, Dn = ace_distribution(trait_col, pollinator_tree)
        pollinator_imputedtraits[!, "L$(coord)_Normal"] = Dn
    end

    pollinator_rec = innerjoin(dropmissing(pollinator_imputedtraits), pollinator_pool; on=:tipNames)#some held-out species may not be in the phylogeny

    #Print out the held-out species that are not in the phylogeny
    missing_in_phylo = setdiff(group_checklist.genus_species, replace.(pollinator_rec.tipNames, "_" => " "))
    if !isempty(missing_in_phylo)
        println("The following genera are not in the phylogeny and will be excluded from predictions: ", join(missing_in_phylo, ", "))
    end
                
    # Sample from phylogeny-inferred distributions
    draws = 20000
    ℒn = Array(pollinator_rec[!, leftnames .* "_Normal"])
    
    # Predict interactions for held-out species
    # Average the latent variable over the 2000 draws
    L_sum = zeros(size(ℒn))
    for i in 1:draws
        L_sum .+= rand.(ℒn)
    end
    L_mean = L_sum ./ draws

    
    #calculate the dot product between L_mean and R
    predicted_dp = L_mean * R

    #Get the observed interaction natrix for the held-out species only
    interaction_matrix_checklist = zeros(Bool, length(pollinator_rec.tipNames), length(plant_in_interaction))
    
    for row in eachrow(group_interaction_df)
        pollinator_idx = findfirst(==(row.genus_species), replace.(pollinator_rec.tipNames, "_" => " "))
        plant_idx = findfirst(==(row.targetTaxonGenusName), plant_in_interaction)
        
        if pollinator_idx !== nothing && plant_idx !== nothing
            interaction_matrix_checklist[pollinator_idx, plant_idx] = true
        end
    end

    #We nned a summary table saving the dot product, predict interaction based on the average threshold, and the observed interaction
    # Get dimensions
    n_pollinators = length(pollinator_rec.tipNames)
    n_plants = length(plant_in_interaction)

    # Create the long-format DataFrame
    summary_df = DataFrame(
        genus_species = repeat(replace.(pollinator_rec.tipNames, "_" => " "), inner=n_plants),
        plant_genus = repeat(plant_in_interaction, outer=n_pollinators),
        dot_product = vec(predicted_dp'),
        predicted_interaction = vec((predicted_dp .>= avg_threshold)'),
        observed_interaction = vec(interaction_matrix_checklist')
    )
    #save the results to a csv file
    CSV.write(
        "Temp/$(group_name)_var_exp_$(v)_prediction_summary.csv",
        summary_df)
end


