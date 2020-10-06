library(shiny)
library(tidyverse)
library(gggenes)
library(RColorBrewer)
library(RCurl)
library(DT)
library(UpSetR)
library(mongolite)
library(shinythemes)
library(shinycssloaders)
library(gggenomes)

# mongod --dbpath /Users/rodtheo/Bioinfo/data-bases/mongo_dbs

#library(svglite)
#library(svgPanZoom)

# p <- ggplot(example_genes, aes(xmin = start, xmax = end, y = molecule, fill = gene)) +
# geom_gene_arrow() +
# facet_wrap(~ molecule, scales = "free", ncol = 1) +
# scale_fill_brewer(palette = getPalette(colourCount))

#x <- getURL("https://raw.githubusercontent.com/rodtheo/entropyDBG/master/multiblast.tsv")

#bgc_table <- read.csv(text=x, stringsAsFactors = FALSE, sep="\t")

db_bgc <- mongo(collection = "bgc", db = "asp_bgc")
db_info <- mongo(collection = "info", db = "asp_bgc")

bgc_table <- read.csv("/Users/rodtheo/Bioinfo/PROJETOS/2019-02-28-NanoporeCompgen/notebook/multiblast_new.tsv", stringsAsFactors = FALSE, sep="\t")

bgc_table %>% mutate(blast_score_num = replace(blast_score, blast_score == 'X', 0)) %>% type_convert() %>% select(blast_score_num) %>% max()


ui <- fluidPage(
    titlePanel("BGC Analysis"),
    #shinythemes::themeSelector(),
    sidebarLayout(
        sidebarPanel(
            radioButtons("refStrain", "What strain is the reference strain ?", db_bgc$distinct({"ref_strain"}), 'NIH'),
            checkboxGroupInput("qryStrains", "Besides the ref strain which strains do you like to exhibit ?",db_bgc$distinct({"ref_strain"}), selected = db_bgc$distinct({"ref_strain"})),
            sliderInput("blastFilter", "Blast Score Filter", value = 1000, min = 0, max = (bgc_table %>% mutate(blast_score_num = replace(blast_score, blast_score == 'X', 0)) %>% type_convert() %>% select(blast_score_num) %>% max())),
            sliderInput("normalizedBlastFilter", "Normalized Blast Score Filter", value = 0.5, min = 0, max = 1.0),
            downloadButton("download_data"),
            actionButton("compare", "Compare BGCs!"),
            radioButtons("source",
                         "Choose to draw",
                         choices = c(
                             "Shared BGCs by ref and selected strains" = "shared",
                             "Unique BGCs in ref" = "unique"
                         )),
            conditionalPanel(
                condition = "input.source == 'shared'",
                uiOutput("typeInput")
            ),
            conditionalPanel(
                condition = "input.source == 'unique'",
                uiOutput("uniqueTypeInput")    
            )
        ),
        mainPanel(
            tabsetPanel(
                tabPanel("Stats", withSpinner(textOutput("textStats"))),
                tabPanel("BGCPlot", withSpinner(plotOutput("coolplot", width = "100%"))),
                #tabPanel("BGCPlot", svgPanZoomOutput(outputId = "coolplot")),
                tabPanel("upset plot", withSpinner(plotOutput("upsetPlot", width = "100%")))
            ),
            textOutput("textN"),
            tabsetPanel(
                tabPanel("BGC reference info:", DT::dataTableOutput("resultsRef")),
                tabPanel("Matched BGCs info:", DT::dataTableOutput("results"))
            )
            
        )
    )
)

server <- function(input, output) {
    
    ## FUNCTIONS
    
    get_core_genes <- function(ref_strain){
        pipeline <- paste0('[
                           {"$match": {"ref_strain": "',ref_strain,'"}},
                           {"$unwind": "$genes"},
                           {"$match": {"genes.description": {"$regex": "biosynthetic .*rule-based-clusters", "$options": "i"}}},
                           {"$project": {"cluster_name": 1, "genes.name": 1, "genes.description": 1, "cluster_type": 1}},
                           {"$group": {"_id": "$cluster_name", "core_genes": {"$addToSet": "$genes.name"}}}
                           ]')
        
        return(db_info$aggregate(pipeline))
    }
    
    check_core_genes <- function(ref_strain, qry_strain, clu_name, matched_clu, table_core_genes, similarity=75) {
        
        
        #print(glimpse(table_core_genes))
        
        ref_cores_set <- table_core_genes %>% dplyr::rename(cluster_name = "_id") %>% dplyr::filter(cluster_name == clu_name)
        
        ref_cores_set <-str_replace_all(ref_cores_set$core_genes[[1]],fixed(".1"),"")
        
        table_similarity <- get_set_of_genes_in_qry_similar(ref_strain, qry_strain, matched_clu)
        
        # print(ref_cores_set)
        # print(table_similarity$homologies$qry_genes$qry_gene)
        
        return( all(is.element(ref_cores_set,  table_similarity$homologies$qry_genes$qry_gene)) )
        
    }
    
    get_set_of_genes_in_qry_similar <- function(ref_strain, qry_strain, matched_clu, sim=75, threshold=0.25) {
        
        pipeline <- paste0('[
                           {"$unwind": "$homologies"},
                           {"$match": {"$and": [{"homologies.score": {"$gte": ',threshold,'}},
                           {"ref_strain": "',ref_strain,'"}, {"homologies.qry_strain": "',qry_strain,'"},                                {"homologies.qry_cluster": "',matched_clu,'"}]}},
                           {"$unwind": "$homologies.qry_genes"},
                           {"$addFields": {"sim_gene": {"$toDecimal": "$homologies.qry_genes.identity"}}},
                           {"$match": {"sim_gene": {"$gte": ',sim,'}}},
                           {"$project": {"_id": 0, "cluster_name": 1, "homologies.qry_cluster": 1, "homologies.qry_genes": 1}}
                           ]')
        
        return(db_bgc$aggregate(pipeline))
    }
    
    
    get_shared_set_two_strains_by_smtype_sim <- function(ref_strain, qry_strain, similarity, show_by_qry_name=FALSE, threshold=0.25){
        
        pipeline <- paste0('[
                           {"$unwind": "$homologies"},
                           {"$match": {"$and": [{"homologies.score": {"$gte": ',threshold,'}}, {"ref_strain": "',ref_strain,'"}, {"homologies.qry_strain": "',qry_strain,'"}]}},
                           {"$group": {"_id": {"clu_name": "$cluster_name"},
                           "set_qry_strains": {"$addToSet": "$homologies.qry_strain"},
                           "set_qry_clusters": {"$push": {"qry": "$homologies.qry_cluster", "score": "$homologies.score"}}}},
                           {"$project": {"score_max": {"$max": "$set_qry_clusters.score"}, "set_qry_clusters": 1}},
                           {"$project": {"max": {"$arrayElemAt": ["$set_qry_clusters", {"$indexOfArray": ["$set_qry_clusters.score",
                           {"$max": "$set_qry_clusters.score"}]}]}}}
                           ]')
        
        all_matches <- db_bgc$aggregate(pipeline) %>% dplyr::rename(ref_clu_name = "_id")
        
        #print(all_matches)
        
        table_core_genes <- get_core_genes(ref_strain)
        ff <- partial(check_core_genes, ref_strain=ref_strain, qry_strain=qry_strain, table_core_genes=table_core_genes, similarity=75)
        zz <- unlist(map2(all_matches$ref_clu_name$clu_name, all_matches$max$qry, ff))
        
        all_matches <- all_matches %>% add_column("match_cores" = zz)
        
        shared_set <- all_matches %>% dplyr::filter(match_cores == TRUE)
        #if( !all(is.element(ref_cores_set,  table_similarity$homologies$qry_genes$qry_gene)) ){
        #    clus_with_cores_not_matched <- c(clus_with_cores_not_matched, clu_name)
        #    matched_with_cores_not_matched <- c(matched_with_cores_not_matched, matched_clu)
        return(shared_set)
    }
    
    get_table_to_draw <- function(ref_clu_name, match_clu_name, table_similarity){
        
        new_table_similarity <- data.frame("cluster_name"=table_similarity$cluster_name, 
                                           "qry_cluster"=table_similarity$homologies$qry_cluster,
                                           "homologous"=table_similarity$homologies$qry_genes$qry_gene,
                                           "name"=str_replace(table_similarity$homologies$qry_genes$qry_gene_original_id, fixed(".1"), ""),
                                           "identity"=table_similarity$homologies$qry_genes$identity,
                                           "coverage"=table_similarity$homologies$qry_genes$coverage,
                                           "evalue"=table_similarity$homologies$qry_genes$evalue)
        
        #x<-c("CH476610.1.region001","contig_7.region002")
        x <- c(ref_clu_name, match_clu_name)
        #pipelines <- map_chr(x, get_cluster_coords)
        tt <- map_df(x, get_cluster_coords)
        
        ttt <- tt %>% mutate(name=str_replace(name, fixed(".1"), ""))
        to_draw_db <- left_join(ttt, new_table_similarity, by="name")
        to_draw_db$homologous <- as.character(to_draw_db$homologous)
        
        to_draw_db <- to_draw_db %>% mutate(sm_genes=case_when(cluster_name.x == ref_clu_name ~ name, TRUE ~ homologous)) %>% mutate(sm_genes=replace_na(sm_genes, "non-homologue"))
        
        return(to_draw_db)
        
    }
    
    get_whole_draw <- function(ref_strain, qry_strain, clu_name, matched_clu){
        table_similarity <- get_set_of_genes_in_qry_similar(ref_strain, qry_strain, matched_clu)
        to_draw_db <- get_table_to_draw(clu_name,matched_clu,table_similarity)
        to_draw_db <- to_draw_db %>% mutate("qry_strain"=rep(qry_strain, length(cluster_name.x)))
        return(to_draw_db)
    }
    
    get_cluster_coords <- function(clu_name){
        pipeline <- paste0('[
                           {"$match": {"cluster_name": "',clu_name,'"}},
                           {"$unwind": "$genes"},
                           {"$project": {"_id": 0, "cluster_name": 1, "name": "$genes.name",
                           "start": "$genes.start", "end": "$genes.end",
                           "strand": "$genes.strand", "description": "$genes.description"}}
                           ]')
        return(db_info$aggregate(pipeline))
    }
    
    plot_using_ggenome <- function(clus_name, g2){
        
        clus_matches <- g2[clus_name][[1]]  %>% distinct(cluster_name.x, name, start, end, .keep_all=TRUE) %>% rename(seq_id=cluster_name.x, feature_id=name) %>% select(seq_id) %>% unique()
        
        clus_matches_str <- str_c(clus_matches$seq_id, collapse=",")[1]
        
        system(paste0("/Users/rodtheo/miniconda2/envs/bioPy3/bin/python /Users/rodtheo/Bioinfo/PROJETOS/2019-02-28-NanoporeCompgen/notebook/script_find_syntenic_regions.py -c ", clus_matches_str))
        
        clus_seqs <- read.csv("/Users/rodtheo/Bioinfo/PROJETOS/2019-02-28-NanoporeCompgen/notebook/seqs_clusters.info.txt", sep="\t", stringsAsFactors = FALSE)
        
        #clus_genes <- 
        
        clus_genes <- g2[clus_name][[1]] %>% distinct(cluster_name.x, name, start, end, .keep_all=TRUE) %>% rename(seq_id=cluster_name.x, feature_id=name) %>% mutate(strand=str_replace_all(strand,"-1","-")) %>% mutate(strand=str_replace_all(strand,"1","+")) %>% select(seq_id, start, end, strand, feature_id)
        
        clus_cogs <-  g2[clus_name][[1]] %>% distinct(cluster_name.x, name, start, end, .keep_all=TRUE) %>% rename(feature_id = name, cluster_id = sm_genes) %>% mutate(cluster_id=replace_na(cluster_id, "non")) %>% select(feature_id, cluster_id)
        
        clus_links <- read_paf("/Users/rodtheo/Bioinfo/PROJETOS/2019-02-28-NanoporeCompgen/notebook/seqs_clusters.paf")
        
        #    clus_genes <- clus_genes %>% mutate(strand=ifelse(seq_id == "contig_8.region007", "+", strand))
        
        p1 <- gggenomes(seqs=clus_seqs, genes=clus_genes, links = clus_links) %>% add_clusters(genes, clus_cogs) + geom_seq() + geom_bin_label() + geom_gene(aes(fill=cluster_id)) + 
            geom_link() 
        return(p1)}
    
    ## SERVER LOGIC
    
    
    filtered <- reactive({
        input$compare
        
        isolate({bgc_table %>% dplyr::filter(ref_strain == input$refStrain, qry_strain %in% input$qryStrains, ref_cluster == input$typeInput) %>% mutate(blast_score_num = replace(blast_score, blast_score == 'X', 0)) %>% group_by(qry_strain) %>% type_convert()})
    })
    
    output$typeInput <- renderUI({
        input$compare
        bgc_sets()
        #qry_gene_tibble <- (filtered() %>% ungroup() %>% select(qry_gene) %>% unique())
        #radioButtons("genes", "What gene ?", qry_gene_tibble$qry_gene)
        #isolate({
            pairwise_clusters <- bgc_sets()$fourth
            #print(head(pairwise_clusters))
            if(input$qryStrains[1] == input$refStrain){
                selected_strain <- input$qryStrains[2]
            } else {
                selected_strain <- input$qryStrains[1]
            }
            
            intersected_df <- pairwise_clusters[[selected_strain]]
            intersected_clus <- intersected_df$ref_clu_name$clu_name
            
            for(i in 1:length(input$qryStrains)){
            
            #for(i in 1:length())
                if(input$qryStrains[i] != input$refStrain){
                    
                    selected_strain <- input$qryStrains[i]
                    selected_clusters_df <- pairwise_clusters[[selected_strain]]
                    selected_clusters <- selected_clusters_df$ref_clu_name$clu_name
                    intersected_clus <- intersect(intersected_clus, selected_clusters)
                }
            }
            
            
            #subset_ref_cluster <- bgc_table %>% filter(ref_strain == input$refStrain) %>% select(ref_cluster) %>% unique()
            selectInput("typeInput", "Ref Cluster Id", choices=intersected_clus)
            # intersected_clus
            # radioButtons("typeInput", "Ref Cluster Id",
            #             choices = sort(subset_ref_cluster$ref_cluster))
        #})
    })
    
    output$uniqueTypeInput <- renderUI({
        input$compare
        
        selectInput("uniqueTypeInput", "Unique BGCs in Ref", choices = sort( bgc_sets()$first ))
    })
    
    
    output$textN <- renderText({
        text_shared <- input$typeInput
        text_uniq <- input$uniqueTypeInput
        
        if(input$source == "shared"){
            print(paste("Shared: ",text_shared))
        }
        else{
            print(paste("Unique: ",text_uniq))
        }
    })
    
    output$textStats <- renderText({
        print("BLEH")
        print("BLAH")
        #cat("This is a test for xxxx\n\nAnd Go on\n\n\ Forever and ever")
    })
    
    rv_synteny <- eventReactive(c(input$compare, input$typeInput), {
        
        db <- bgc_sets()$third
        
        colourCount <- length(unique(db[[which(names(db) == input$typeInput)]]$homologous))
        getPalette <- colorRampPalette(brewer.pal(12, "Set3"))
        
        
        plot_using_ggenome(input$typeInput, db)
        
       # subset_data <- db[[which(names(db) == input$typeInput)]] %>% distinct(cluster_name.x, name, start, end, .keep_all=TRUE) %>% dplyr::filter(qry_strain %in% input$qryStrains)
        
        
        
        
        #+ scale_fill_manual(values = c(colorRampPalette(brewer.pal(12, "Set3"))(colourCount-1), "#FFFFFF")) + theme_genes() + ylab("")

    })
    
    output$coolplot <- renderPlot({
        
        rv_synteny()
        
        #subset_data <- rv_synteny()
        
        
        
        
        #ggplot(as_tibble(subset_data), aes(xmin = start, xmax =end, y = cluster_name.x, forward=strand, fill=sm_genes, label=homologous)) + geom_gene_arrow(arrowhead_height = unit(3, "mm"), arrowhead_width = unit(1, "mm")) +  geom_gene_label() + facet_wrap(~ cluster_name.x, scales = "free", ncol = 1) + theme_genes()
    
        # 
        # ref_strain <- "NIH"
        # df_ref <- filtered() %>% filter(blast_score_num == max(blast_score_num))
        # 
        # df_ref_to_plot <- df_ref %>% select(qry_gene, qry_gene_start, qry_gene_end, qry_gene_strand, qry_strain)
        # 
        # # df_ref_to_plot$qry_strain <- as.factor(df_ref_to_plot$qry_strain)
        # # qry_strain_levels <- levels(df_ref_to_plot$qry_strain)
        # # 
        # # # remove the ref_strain from levels vector and put it in front of updated level
        # # update_qry_strain_levels <- c(ref_strain, qry_strain_levels[!(qry_strain_levels == ref_strain)])
        # # 
        # # levels(df_ref_to_plot$qry_strain) <- update_qry_strain_levels
        # 
        # df_ref_to_plot$qry_strain <- fct_relevel(as.factor(df_ref_to_plot$qry_strain), "NIH")
        # 
        # colourCount <- length(unique(df_ref_to_plot$qry_gene))
        # getPalette <- colorRampPalette(brewer.pal(12, "Set3"))
        # 
        # ggplot(df_ref_to_plot, aes(xmin = qry_gene_start, xmax = qry_gene_end, y = qry_strain, fill = qry_gene, forward=qry_gene_strand, label=qry_gene)) + geom_gene_arrow(arrowhead_height = unit(3, "mm"), arrowhead_width = unit(1, "mm")) +  geom_gene_label() + facet_wrap(~ qry_strain, scales = "free", ncol = 1) + scale_fill_manual(values = c(colorRampPalette(brewer.pal(12, "Set3"))(colourCount-1), "#FFFFFF")) + theme_genes() + ylab("")
    })
    
    output$resultsRef <- DT::renderDataTable({
        df_ref_ref <- filtered() %>% dplyr::filter(ref_strain == input$refStrain, qry_strain == input$refStrain) %>% select(qry_gene, qry_gene_start, qry_gene_end, qry_gene_strand, qry_strain, gene_function)
        #datatable(df_ref, rownames = FALSE, filter="top", options = list(pageLength = 5, scrollX=T) )
        df_ref_ref
    }, extensions = 'Buttons',
    options = list(dom = 'Bfrtip',
                   buttons = c('copy', 'csv', 'excel', 'pdf', 'print')))
    
    output$download_data <- downloadHandler(
        filename = "bgc_table.csv",
        content = function(file) {
            data <- filtered() %>% dplyr::filter(ref_strain == input$refStrain, qry_strain == input$refStrain) %>% select(qry_gene, qry_gene_start, qry_gene_end, qry_gene_strand, qry_strain, gene_function)
            write.csv(data, file, row.names = FALSE)
        }
    )
    
    output$results <- DT::renderDataTable({
        df_ref <- filtered() %>% dplyr::filter(blast_score_num >= input$blastFilter, normalized_score >= input$normalizedBlastFilter) %>% select(qry_strain, qry_gene, coverage, identity, evalue, qry_gene_start, qry_gene_end, qry_gene_strand,  blast_score_num, normalized_score, gene_function, qry_cluster)
        #datatable(df_ref, rownames = FALSE, filter="top", options = list(pageLength = 5, scrollX=T) )
        df_ref
    }, extensions = 'Buttons',
    fillContainer = FALSE,
    options = list(dom = 'Bfrtip',
                   buttons = c('copy', 'csv', 'excel', 'pdf', 'print')))
    
    
    bgc_sets <- reactive({
        input$compare
        isolate({
        ref_strain <- input$refStrain
        qry_strains <- input$qryStrains
        #qry_strains <- c("45A", "IFO6365")
        
        qry_strain <- "antismash_ATCC20542_gbk_reano"
        
        z <- get_shared_set_two_strains_by_smtype_sim(ref_strain, qry_strain, 85, threshold = 0.25)
        
        shared_set <- z$ref_clu_name$clu_name
        
        results_list <- list()
        results_list$first <- shared_set
        results_list$second <- shared_set
        
        
        
        pairwise_sets <- list()
        to_draw_ncbi <- c()
        strains_order <- c()
        for(i in 1:length(qry_strains)){
            if (qry_strains[i] != input$refStrain){
            #if (qry_strains[i] != ref_strain){
            print(qry_strains[i])
            z <- get_shared_set_two_strains_by_smtype_sim(ref_strain, qry_strains[i], 75, threshold = 0.25)
            pairwise_sets[[qry_strains[i]]] <- z
            strains_order <- c(strains_order, qry_strains[i])
            
            partial_get_whole_draw <- partial(get_whole_draw, ref_strain=ref_strain, qry_strain=qry_strains[i])
            
            to_draw_ncbi_qry <- map2(z$ref_clu_name$clu_name, z$max$qry, partial_get_whole_draw)
            names(to_draw_ncbi_qry) <- z$ref_clu_name$clu_name
            
            to_draw_ncbi <- c(to_draw_ncbi, to_draw_ncbi_qry)
            }}
        
        
        glist <- to_draw_ncbi
        g2 <- map(split(glist, names(glist)), bind_rows)
        
        results_list$third <- g2
        
        results_list$fourth <- pairwise_sets
        
        #print(head(results_list$third))
        
        return(results_list)})
        # print(listInput)
        # print(names(listInput))
    })
    
    

    
    output$upsetPlot <- renderPlot({
        
        # print(names(listInput))
        #names(listInput) <- group_keys(y)$qry_strain
        #print(listInput)
        #upset(fromList(listInput), sets.bar.color = "#56B4E9", order.by = "freq", nsets = length(input$qryStrains))
        p <- upset(fromExpression(bgc_sets()$second), sets.bar.color = "#56B4E9", order.by = "freq", nsets = length(input$qryStrains))
        p
    })
}

shinyApp(ui, server)

