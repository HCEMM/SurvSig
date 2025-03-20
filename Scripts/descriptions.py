import streamlit as st


def rousseaux_surv_type():
    text = ':blue[Plots Need ] :red[Refresh if You Change the Survival Type]'

    return text


def site_description():
    logo_col = st.columns(5)
    logo_col2 = st.columns(3)
    cont_login = st.container(border=True)
    logo_col[2].image("style_items/coming.svg", use_container_width=True)
    logo_col2[1].video("videos/WelcometoS.mp4", autoplay=True, muted=True, loop=True)
    cont_login.header(":blue[Description]")
    cont_login.write("Welcome to our advanced cancer data analysis tool! Powered by machine learning algorithms, "
                     "our application specializes in clustering gene expression data, revealing intricate patterns "
                     "and insights into various cancer types. Beyond clustering, the tool also provides "
                     "comprehensive statistical analyses to give a deeper understanding of the data. Through the "
                     "synergy of genomics and cutting-edge technology, we aim to advance the landscape of cancer "
                     "research.")
    #cont_login.header(":blue[Special Thanks]")
    #cont_login.write("We deeply appreciate the time and effort our volunteer testers dedicated to refining this "
    #                 "application. Your feedback and dedication have been instrumental in bringing this tool to "
    #                "the forefront. Together, we are pioneering the future of cancer research. Thank you!")
    cont_login.warning(":red[Note: ] :blue[If you ] :red[refresh ] :blue[the webpage, the user will be]"
                       " :red[logged out]!")


def mds_adv():
    text = (":blue[Multidimensional scaling. (MDS).]\n "
            "\n:red[-Use Metric MDS:] If True, perform metric MDS; otherwise, perform nonmetric MDS. When False "
            "(i.e. non-metric MDS), dissimilarities with 0 are considered as missing values.\n"
            "\n:red[-Maximum Number of Iterations:] Maximum number of iterations of the SMACOF algorithm for a single run.\n"
            "\n:red[-Epsilon Value:] Relative tolerance with respect to stress at which to declare convergence.\n")

    return text


def nmf_adv():
    text = (":blue[Non-Negative Matrix Factorization (NMF).]\n "
            "Find two non-negative matrices, i.e. matrices with all non-negative elements, (W, H) whose product "
            "approximates the non-negative matrix X. This factorization can be used for example for "
            "dimensionality reduction, source separation or topic extraction.\n"
            "\n:red[-Maximum Number of Iterations: ] Maximum number of iterations before timing out.\n"
            "\n:red[-Tolerance of Stopping Condition: ] Tolerance of the stopping condition.\n"
            "\n:red[-Method of Initialize: ] Method used to initialize the procedure.\n"
            "\n:red[-Beta Loss: ] Beta divergence to be minimized, measuring the distance between X and the dot "
            "product WH. Note that values different from ‘frobenius’ (or 2) and ‘kullback-leibler’ (or 1) lead to "
            "significantly slower fits.\n"
            "\n:red[-Select Numerical Solver: ] cd’ is a Coordinate Descent solver and ‘mu’ is a Multiplicative Update "
            "solver. \n"
            "\n:red[-Randomize the Order of Coordinates: ] If true, randomize the order of coordinates in the CD "
            "solver.\n")

    return text


def pca_adv():
    text = (":blue[Principal component analysis (PCA)]\n "
            "\n:red[-Relative Tolerance: ] Tolerance for singular values computed by svd_solver == ‘arpack’. "
            "Must be of range [0.0, infinity).\n"
            "\n:red[-Use Whiten: ] When True (False by default) Whitening will remove some information from the "
            "transformed signal (the relative variance scales of the components) but can sometime improve the "
            "predictive accuracy of the downstream estimators by making their data respect some hard-wired "
            "assumptions. \n")

    return text


def umap_adv():
    text = (":blue[Uniform Manifold Approximation and Projection for Dimension Reduction (UMAP)]\n "
            "\n Uniform Manifold Approximation and Projection (UMAP) is a dimension reduction technique that can be "
            "used for visualisation similarly to t-SNE, but also for general non-linear dimension reduction. \n"
            "\n UMAP is a fairly flexible non-linear dimension reduction algorithm. It seeks to learn the manifold "
            "structure of your data and find a low dimensional embedding that preserves the essential topological "
            "structure of that manifold. \n"
            "\n:red[-Minimal Distance Between Points: ] The min_dist parameter controls how tightly UMAP is allowed "
            "to pack points together. It, quite literally, provides the minimum distance apart that points are allowed "
            "to be in the low dimensional representation. This means that low values of min_dist will result in "
            "clumpier embeddings. This can be useful if you are interested in clustering, or in finer topological "
            "structure. Larger values of min_dist will prevent UMAP from packing points together and will focus on the "
            "preservation of the broad topological structure instead.\n"
            "\n:red[-Number of Neighbours: ] This parameter controls how UMAP balances local versus global structure "
            "in the data. It does this by constraining the size of the local neighborhood UMAP will look at when "
            "attempting to learn the manifold structure of the data. This means that low values of n_neighbors will "
            "force UMAP to concentrate on very local structure (potentially to the detriment of the big picture), "
            "while large values will push UMAP to look at larger neighborhoods of each point when estimating the"
            " manifold structure of the data, losing fine detail structure for the sake of getting the broader of "
            "the data.\n"
            "\n:red[-Select Metric: ] This controls how distance is computed in the ambient space of the input data.\n")

    return text


def phate_adv():
    text = (":blue[Potential of Heat-diffusion for Affinity-based Trajectory Embedding (PHATE)]\n "
            "\n PHATE is a tool for visualizing high dimensional data. PHATE uses a novel conceptual framework for "
            "learning and visualizing the manifold to preserve both local and global distances. \n"
            "\n:red[-Number of Neighbours: ] Number of nearest neighbors on which to build kernel\n"
            "\n:red[-Number of Principal Components: ] Number of principal components to use for calculating "
            "neighborhoods. For extremely large datasets, using n_pca < 20 allows neighborhoods to be calculated in "
            "roughly log(n_samples) time.\n")

    return text


def tsne_adv():
    text = (":blue[T-distributed Stochastic Neighbor Embedding (t-SNE)]\n "
            "\n t-SNE is a tool to visualize high-dimensional data. It converts similarities between data points to "
            "joint probabilities and tries to minimize the Kullback-Leibler divergence between the joint probabilities "
            "of the low-dimensional embedding and the high-dimensional data. t-SNE has a cost function that is "
            "not convex, i.e. with different initializations we can get different results. \n"
            "\n:red[-Perplexity of t-SNE: ] The perplexity is related to the number of nearest neighbors that is used "
            "in other manifold learning algorithms. Larger datasets usually require a larger perplexity. "
            "Consider selecting a value between 5 and 50. Different values can result in significantly "
            "different results. The perplexity must be less than the number of samples. \n"
            "\n:red[-Learning Rate: ] The learning rate for t-SNE is usually in the range [10.0, 1000.0]. "
            "If the learning rate is too high, the data may look like a ‘ball’ with any point approximately "
            "equidistant from its nearest neighbours. If the learning rate is too low, most points may look "
            "compressed in a dense cloud with few outliers. If the cost function gets stuck in a bad local minimum "
            "increasing the learning rate may help.\n"
            "\n:red[-Distance Metric: ] The metric to use when calculating distance between instances in a "
            "feature array.\n"
            "\n:red[-Number of Iterations: ] Maximum number of iterations for the optimization. "
            "Should be at least 250. \n"
            "\n:red[-Initialisation of Embedding Method: ] Initialization of embedding. PCA initialization cannot be"
            " used with precomputed distances and is usually more globally stable than random initialization.\n"
            )

    return text


def nmf_clust_adv():
    text = (":blue[Non-Negative Matrix Factorization Clustering (NMF Clustering)]\n "
            "\n NMF Clustering factorizes a non-negative input matrix into non-negative factors. The algorithm has an "
            "inherent clustering property and has been gaining attention in various fields especially in "
            "biological data analysis. \n"
            "\n:red[-[Number of Iterations]: ]  Determine how many iterations we want to run.\n"
            "\n:red[-Number of Trials: ]  Determine how many trials we want to run.\n"
            "\n:red[-Lamb's Value: ]  Hyper-parameter for the Integrative NMF algorithm that controls the rate of "
            "learning\n"
            )

    return text


def corr_adv():
    text = (":blue[Correlation]\n "
            "\n Pearson correlation coefficient and p-value for testing non-correlation. The Pearson correlation "
            "coefficient measures the linear relationship between two datasets. Calculate a Spearman correlation "
            "coefficient with associated p-value. The Spearman rank-order correlation coefficient is a nonparametric "
            "measure of the monotonicity of the relationship between two datasets.\n"
            "\n:red[-Correlation Method: ] You can use Pearson for parametric data, and you can use Spearman for "
            "nonparametric data  \n"
            )

    return text


def gene_list_warning():
    return st.warning(
           ":blue[Please ] :red[upload your gene list] :blue[or ] :red[you can make one]! "
           ":blue[Check the ] :red[Gene List TAB ] :blue[or ] :red[Gene Set Finder TAB] "
           ":blue[for more information!]")


def demo_gene_list_help():
    text = (":blue[Demo Gene List for SurvSig]\n"
            "\nThe *NE50 gene list* serves as a demonstration dataset within the **SurvSig** application. This curated list includes "
            "**50 genes**: 25 associated with **neuroendocrine (NE)** characteristics and 25 linked to **non-neuroendocrine (non-NE)** traits. "
            "Designed as a general example, the NE50 gene list provides a foundation for exploring gene expression and survival analysis in "
            "**Neuroendocrine Lung Tumors**, including **Small Cell Lung Cancer (SCLC)**.\n"
            "\n:red[This dataset highlights] key gene expression patterns and supports users in understanding the functionality of SurvSig, "
            "offering insights into NE and non-NE classifications and their potential implications for clinical research and tumor biology.\n"
            )
    return text



def opt_adv():
    text = (":blue[Ordering Points To Identify the Clustering Structure (OPTICS)]\n "
            "\n OPTICS, closely related to DBSCAN, finds core sample of high density and expands clusters from them. "
            "Unlike DBSCAN, keeps cluster hierarchy for a variable neighborhood radius. \n"
            "\n:red[-Minimum Samples: ] The number of samples in a neighborhood for a point to be considered as a "
            "core point. Also, up and down steep regions can’t have more than min_samples consecutive non-steep "
            "points. \n"
            "\n:red[-Distance Metric: ] Metric to use for distance computation \n"
            "\n:red[-Initialisation Method: ] The extraction method used to extract clusters using the calculated "
            "reachability and ordering. \n"
            "\n:red[-Cluster Selection Epsilon Value: ] The maximum distance between two samples for one to be "
            "considered as in the neighborhood of the other. By default it assumes the same value as max_eps. "
            "Used only when cluster_method='dbscan'. \n"
            )

    return text


def hdb_adv():
    text = (":blue[Hierarchical Density-Based Spatial Clustering of Applications with Noise (HDBSCAN)]\n "
            "\n Cluster data using hierarchical density-based clustering. Performs DBSCAN over varying epsilon values "
            "and integrates the result to find a clustering that gives the best stability over epsilon. This allows"
            " HDBSCAN to find clusters of varying densities (unlike DBSCAN), and be more robust to parameter "
            "selection. \n"
            "\n:red[-Alpha Value: ] A distance scaling parameter as used in robust single linkage. \n"
            "\n:red[-Metric Method: ] The metric to use when calculating distance between instances in a "
            "feature array. \n"
            "\n:red[-Minimum Number of Samples in Each Cluster: The minimum number of samples in a group for that "
            "group to be considered a cluster; groupings smaller than this size will be left as noise.]  \n"
            "\n:red[-Cluster Selection Epsilon Value: ] A distance threshold. Clusters below this value will be "
            "merged. \n"
            "\n:red[-Select the Cluster Selection Method: ] The method used to select clusters from the condensed tree. "
            "The standard approach for HDBSCAN* is to use an Excess of Mass ('eom') algorithm to find the most "
            "persistent clusters. Alternatively you can instead select the clusters at the leaves of the "
            "tree – this provides the most fine grained and homogeneous clusters. \n"
            )

    return text


def dtc_adv():
    text = (":blue[dynamicTreeCut]\n "
            "\n A variable height branch pruning technique for dendrograms produced by hierarchical clustering. \n"
            "\n:red[-Distance Metric: ] Pairwise distances between observations in n-dimensional space. \n"
            "\n:red[-Linkage Method: ] The linkage algorithm to use. Perform hierarchical/agglomerative clustering.  \n"
            "\n:red[-Minimum Cluster Size: ] Each cluster may contain how many samples. \n"
            )

    return text


def som_adv():
    text = (":blue[Self Organizing Maps (SOM)]\n "
            "\n  SOM is a type of Artificial Neural Network able to convert complex, nonlinear statistical "
            "relationships between high-dimensional data items into simple geometric relationships on a low-dimensional "
            "display. \n"
            "\n:red[-Neighbours Function Type: ]  Function that weights the neighborhood of a position in the map. \n"
            "\n:red[-Sigma Value: ]  Spread of the neighborhood function, needs to be adequate to the "
            "dimensions of the map. \n"
            "\n:red[-Learning Rate: ]  Initial learning rate. \n"
            "\n:red[-Number of Epochs: ]  Number of epochs. An epoch in machine learning means one complete pass of "
            "the training dataset through the algorithm.\n"
            "\n:red[-Choose the Activation Distance: ]  Distance used to activate the map. \n"
            "\n:red[-Select the Topology of the Map: ] Topology of the map. \n"
            )

    return text


def kmeans_adv():
    text = (":blue[K-Means Clustering]\n "
            "\n K-means is an unsupervised learning method for clustering data points. The algorithm iteratively "
            "divides data points into K clusters by minimizing the variance in each cluster. \n"
            "\n:red[-K-means Initialisation Method: ]  Method for initialization:\n"
            "\n -‘k-means++’ : selects initial cluster centroids using sampling based on an empirical probability "
            "distribution of the points’ contribution to the overall inertia. This technique speeds up convergence. "
            "The algorithm implemented is “greedy k-means++”. It differs from the vanilla k-means++ by making several"
            " trials at each sampling step and choosing the best centroid among them. \n"
            "\n ‘random’: choose n_clusters observations (rows) at random from data for the initial centroids. \n"
            "\n:red[-Relative Tolerance: ]  Relative tolerance with regards to Frobenius norm of the difference in "
            "the cluster centers of two consecutive iterations to declare convergence.\n"
            "\n:red[-Max iteration: ] Maximum number of iterations of the k-means algorithm for a single run. \n"
            "\n:red[-K-means Algorithm to Use: ]  K-means algorithm to use. The classical EM-style algorithm "
            "is 'lloyd'. The 'elkan' variation can be more efficient on some datasets with well-defined clusters,"
            " by using the triangle inequality. However it’s more memory intensive due to the allocation of an extra"
            " array of shape (n_samples, n_clusters).\n"
            )

    return text


def gmm_adv():
    text = (":blue[Gaussian Mixture Model (GMM)]\n "
            "\n Representation of a Gaussian mixture model probability distribution. This class allows to estimate"
            " the parameters of a Gaussian mixture distribution. \n"
            "\n:red[-Initialisation Methods: ]  The method used to initialize the weights, the means and the "
            "precisions. String must be one of:\n"
            "\n -‘kmeans’ : responsibilities are initialized using kmeans. \n"
            "\n -‘k-means++’ : use the k-means++ method to initialize. \n"
            "\n -‘random’ : responsibilities are initialized randomly. \n"
            "\n -‘random_from_data’ : initial means are randomly selected data points. \n"
            "\n:red[-Convariance Type: ]  String describing the type of covariance parameters to use. Must be one of:\n"
            "\n -‘full’: each component has its own general covariance matrix. \n"
            "\n -‘tied’: all components share the same general covariance matrix. \n"
            "\n -‘diag’: each component has its own diagonal covariance matrix. \n"
            "\n -‘spherical’: each component has its own single variance. \n"
            "\n:red[-Convergence Threshold: ]  The convergence threshold. EM iterations will stop when the lower "
            "bound average gain is below this threshold.\n"
            "\n:red[-Non-negative Regularization: ]  Non-negative regularization added to the diagonal of covariance."
            " Allows to assure that the covariance matrices are all positive.\n"
            )

    return text


def agglo_adv():
    text = (":blue[Agglomerative Clustering]\n "
            "\n Recursively merges pair of clusters of sample data; uses linkage distance. \n"
            "\n:red[-Metric of Agglomerative Clustering: ]  Metric used to compute the linkage. Can be “euclidean”, "
            "“l1”, “l2”, “manhattan”, “cosine”, or “precomputed”. If linkage is “ward”, only “euclidean” is accepted. "
            "If “precomputed”, a distance matrix is needed as input for the fit method.\n"
            "\n:red[-Linkage Method of Agglomerative Clustering: ] Which linkage criterion to use. "
            "The linkage criterion determines which distance to use between sets of observation. The algorithm "
            "will merge the pairs of cluster that minimize this criterion. \n"
            "\n -‘ward’ minimizes the variance of the clusters being merged. \n"
            "\n -‘average’ uses the average of the distances of each observation of the two sets. \n"
            "\n -‘complete’ or ‘maximum’ linkage uses the maximum distances between all observations of the two sets. \n"
            "\n -‘single’ uses the minimum of the distances between all observations of the two sets. \n"
            )

    return text


def tcga_description():
    text = ("\n :red[The Cancer Genome Atlas (TCGA) Project]\n"
            "\n :blue[Overview:]\n"
            "\nTCGA is a comprehensive cancer genomics dataset that encompasses data from over 11,000 patients across 33 different tumor types. "
            "This project, launched by the National Institutes of Health (NIH), has significantly contributed to our understanding of cancer by providing "
            "detailed molecular profiles, including gene expression, mutations, copy number variations, and epigenetic modifications.\n"
            "\n :blue[Tumor Types:]\n"
            "\nTCGA includes a wide range of cancers, from common types like :red[BRCA (Breast Invasive Carcinoma)], :red[LUAD (Lung Adenocarcinoma)], "
            "and :red[COAD (Colon Adenocarcinoma)] to rarer forms such as :red[ACC (Adrenocortical Carcinoma)], :red[UVM (Uveal Melanoma)], and "
            ":red[SARC (Sarcoma)]. Each tumor type is represented by hundreds of samples, providing a robust dataset for comprehensive analyses.\n"
            "\n :blue[Data Available:]\n"
            "\nThe TCGA dataset offers multi-dimensional data, including transcriptomics, genomics, epigenomics, and proteomics. This data has been instrumental "
            "in identifying novel biomarkers, understanding the molecular basis of cancer, and advancing personalized medicine approaches.\n"
            "\n :blue[Publications and Impact:]\n"
            "\nThe TCGA project has led to numerous high-impact publications and continues to be a vital resource for cancer research worldwide. "
            "One of the seminal publications from TCGA is available at [DOI: 10.1038/ng.2764](https://doi.org/10.1038/ng.2764), detailing the comprehensive "
            "analysis of cancer genomes.\n"
            "\n :blue[* Data Preprocessing Information *]"
            "\n :red[* Only samples with transcriptomic data were included in the analysis. For survival analysis, only samples with both transcriptomic and survival data were used. *]\n"
            "\n :red[* Gene names were standardized to HUGO symbols if the original study used a different nomenclature. Genes that could not be precisely mapped to HUGO symbols were removed from the dataset. *]\n"
            "\n :red[* For more detailed information on the specific cancer types included in TCGA, you can visit the official resource page: *]\n"
            "[TCGA Study Abbreviations](https://gdc.cancer.gov/resources-tcga-users/tcga-code-tables/tcga-study-abbreviations)\n")

    return text


def exp_type():
    text = ('First, select the type of gene expression database you would like to explore. You have two main options: '
            ':red["Neuroendocrine Lung Tumor"] and :red["TCGA"].\n\n'
            ':red["Neuroendocrine Lung Tumor"] includes seven distinct cohorts specifically focused on various types of neuroendocrine lung tumors, such as '
            ':blue[George-SCLC], :blue[Anish-SCLC], :blue[Jiang-SCLC], :blue[Liu-SCLC], :blue[George-LCNEC], :blue[Fernandez-Carcinoid], '
            ':blue[Alcala-Carcinoid] and :blue[Rousseax-Mixed]. These cohorts provide valuable insights into the gene expression profiles of different subtypes of small cell lung cancer '
            'and related neuroendocrine tumors.\n\n'
            'On the other hand, :red["TCGA"] (The Cancer Genome Atlas) encompasses a much broader range of 33 tumor types, offering gene expression data '
            'across a wide spectrum of cancers, including common types such as breast, lung, and colon cancer, as well as rarer forms. This comprehensive '
            'dataset allows for cross-cancer comparisons and a deeper understanding of gene expression patterns in various malignancies.')

    return text


def download():
    text = (':blue[Add name for your ] :red[results!]\n')

    return text



def gl_text():
    text = ('The :red[mapping of gene names ] is necessary in order to get consistent names in our work. You can '
            ':red[select a suggested gene name or delete ] it from the list.')

    return text


def sub_text():
    text = ('The Subtype Select feature allows you to filter data using :red["OR" ] logic for inclusive results across'
            ' any of your chosen subtypes or :red["AND" ] logic to exclusively display entries that meet all selected '
            'criteria. Use :red["OR" ] to broaden your search and :red["AND" ] to narrow it down to specific entries. '
            'By toggling the logic options, you can easily switch between a wider dataset and a more focused one. '
            'Check the desired subtypes, '
            'set your logic, and the system will display the corresponding results.')

    return text


def ssgsea_surv_intro():
    text = ('Select the checkbox to :red[start ssGSEA Survival analysis. ]'
            "ssGSEA, short for single-sample :blue[Gene Set Enrichment Analysis], "
            "quantifies gene set activity in individual "
            "samples. By sorting samples into groups based on ssGSEA scores and analyzing their survival data, "
            "researchers can identify gene expression's impact on patient outcomes, offering "
            "insights into disease prognosis.")

    return text


def sep_method():
    text = ('\n:blue[Automatic: ] The app identifies the most statistically significant cutoff for separating samples, '
            'streamlining the process for the user. \n'
            '\n:blue[Median: ] This option splits the sample set into two groups based on the median ssGSEA score. \n'
            '\n:blue[Percentage Decomposition: ] Samples are divided into two groups at a specific percentile cutoff '
            'defined by the user.\n '
            '\n:blue[Percentage Decomposition ] (lower and upper limit): This creates three groups of samples by '
            'specifying both lower and upper percentile cutoffs, enabling more detailed stratification. \n'
            '\n:blue[ssGSEA Values Cutoff: ] Users can manually enter a specific ssGSEA score to separate '
            'the samples into two groups.')

    return text


def sep_method_sg():
    text = ('\n:blue[Automatic: ] The app identifies the most statistically significant cutoff for separating samples, '
            'streamlining the process for the user. \n'
            '\n:blue[Median: ] This option splits the sample set into two groups based on the median ssGSEA score. \n'
            '\n:blue[Percentage Decomposition: ] Samples are divided into two groups at a specific percentile cutoff '
            'defined by the user.\n '
            '\n:blue[Percentage Decomposition ] (lower and upper limit): This creates three groups of samples by '
            'specifying both lower and upper percentile cutoffs, enabling more detailed stratification. \n'
            '\n:blue[Expression Values Cutoff: ] Users can manually enter a specific expression level to separate '
            'the samples into two groups.')

    return text


def perc_int():
    text = ('This controls the span of sample distribution to be analyzed, allowing selection of the most relevant'
            ' sections based on ssGSEA score percentiles.')

    return text


def sg_intro():
    text = ("The :blue[Single Gene Survival ] tab lets you see the influence of "
            ":red[one gene's expression on survival.] "
            " Just select a gene and the app provides an optimal expression threshold, along with"
            " visual plots for expression levels, survival comparison, and a multivariate analysis"
            " for a concise yet comprehensive survival impact assessment.")

    return text


def go_intro():
    text = ('The :blue[Gene Ontology ] feature in this tool allows for the functional analysis of genes.'
            ' You can select a :red[set of genes or gene clusters ] to analyze their roles and relationships. ')

    return text


def go_geneset():
    text = ('For the :blue["Gene Set" ] option, choose :Red["Gene List" ] to analyze a custom selection of genes, '
            'or :red["Cluster" ] to study a pre-defined group of genes, allowing '
            'for tailored or broader functional insights.')

    return text


def multi_intro():
    text = (":blue[Multivariate ] analysis helps to see how different factors work together to impact a result. "
            "The :blue[Chi-square test ] checks if there's a clear link between two things we can count or categorize.")

    return text


def dataset_ana_intro():
    text = (
        ':blue[Gene Set Finder]: This feature leverages a variety of :red[machine learning] and :red[statistical methods] '
        'to help you identify significant genes within a selected cohort. By analyzing gene expression data, it identifies key genes that '
        ':red[characterize] specific aspects of your dataset, such as cancer subtypes, stages, or other clinical outcomes. This tool is designed to '
        'discover :red[biomarkers] and :red[predictive gene sets] that can be used for diagnosis, prognosis, and personalized treatment strategies. '
        'The integration of clinical data further enhances the precision of the analysis, making it a powerful resource for cancer research and clinical applications.')
    return text


def go_fun():
    text = ('This section allows you to use several :blue[predefined statistical methods] '
            'for :red[enrichment analysis. ] '
            'The choice of function typically depends on the :red[type of enrichment analysis you want to conduct. ] \n'
            '\n":blue[enrichGO": ] To conduct GO enrichment analysis using over-representation analysis.\n'
            '\n":blue[enrichKEGG": ] For KEGG pathway enrichment analysis.\n')

    return text


def p_value_cutoff():
    text = ':red[Adjusted P value ] cutoff on :blue[enrichment tests ] to report'

    return text


def q_value_cutoff():
    text = ':red[Q value cutoff ] on :blue[enrichment tests ] to report as :red[significant.]'

    return text


def corr_meth_go():
    text = ':red[One ] of :blue["holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"]'

    return text


def max_gene():
    text = ':red[Maximal size ] of genes :blue[annotated ] for testing'

    return text


# Function for Description argument
def gene_ontology_title():
    text = ':red[Title]: :blue[Title] of the :red[Plot].'
    return text

# Function for Max Category argument
def max_category_description():
    text = ':red[Max Category]: Specifies the :blue[maximum number] of categories displayed in the analysis results.'
    return text

# Function for OrderBy argument
def order_by_description():
    text = ':red[Order By]: Determines the :blue[criterion] for sorting the pathways, such as by :blue[gene ratio] or :blue[p-value].'
    return text


def sub_ont():
    text = (
        '\n:blue[BP (Biological Process): ] This ontology refers to biological objectives to which the gene or gene '
        'product contributes. For example, biological processes include '
        'photosynthesis, cell cycle, or immune response.\n'
        '\n:blue[CC (Cellular Component): ] This ontology describes the '
        'location within the cell where a gene product '
        'is active, including structures such as the nucleus, cytoplasm, or membrane.\n'
        '\n:blue[MF (Molecular Function): ] This ontology covers the biochemical activities at the molecular '
        'level performed by gene products, such as enzyme activities, binding functions, or receptor activities.\n')

    return text


def act_btn():
    text = 'Click for :red[accept ] the settings!'

    return text


def sep_step():
    text = ('This adjusts the increment between the percentile markers, fine-tuning the range within which the '
            'Spearman correlation is calculated for sample separation.')

    return text


def ssgsea_comp_intro():
    text = ('The :blue[ssGSEA ] function in this tab helps you understand gene activity within '
            'your samples by :red[comparing enrichment scores. ] You can assess the relationship '
            'between :red[different gene sets or between a gene set and a single gene, ] revealing '
            'their correlation in an intuitive way. Clustering genes is recommended.')

    return text


def gene_comp_into():
    text = (':blue[Gene Compare] section lets you :red[compare the expression of two genes ] '
            'to see if they change together across different '
            'samples, giving you insights into their potential interaction or common pathways in your study.')

    return text


def describe_feat_selector_meth2():
    text = (
        ':blue[Feature Selector Method]: Choose the method for :red[selecting features] from the data to identify the most '
        'informative genes that can predict cancer outcomes such as NE, EMT, histology, grade, etc.\n'
        '\n:red[Options:]\n'
        '- :blue[PCA (Principal Component Analysis)]: A technique that reduces the dimensionality of the data while retaining the most '
        'variance, helping in identifying the most significant genes.\n'
        '- :blue[ICA (Independent Component Analysis)]: Separates a multivariate signal into additive, independent components, useful for '
        'identifying independent sources of variation in gene expression.\n'
        '- :blue[FA (Factor Analysis)]: Identifies underlying relationships between variables by modeling observed variables and their '
        'correlations, useful for clustering genes with similar expression patterns.\n'
        '- :blue[SVD (Singular Value Decomposition)]: A matrix factorization method that decomposes the data into singular vectors and values, '
        'useful for noise reduction and identifying key gene expression patterns.\n'
        '- :blue[Standard Deviation]: Selects genes based on their standard deviation, prioritizing genes with the most variation across samples.\n'
        '- :blue[NMF (Non-Negative Matrix Factorization)]: Factorizes the data matrix into non-negative matrices, useful for parts-based '
        'representation of gene expression data and uncovering hidden structures.')
    return text


def describe_rank_method2():
    text = (
        ':blue[Ranking Method]: Select the method to :red[rank the features] based on their importance in predicting cancer outcomes.\n'
        '\n:red[Options:]\n'
        '- :blue[Sum Square]: Ranks genes based on the sum of squared loadings across components. Use this when you want to emphasize genes '
        'that contribute most to multiple principal components, ensuring a holistic view of gene importance.\n'
        '- :blue[Average Absolute]: Ranks genes based on the average absolute loadings across components. This method is useful when you want '
        'to identify genes that have a consistent influence across different components, providing a balanced selection.\n'
        '- :blue[Max Absolute]: Ranks genes based on the maximum absolute loading in any component. Choose this method to highlight genes '
        'that have a strong influence in at least one component, ideal for identifying highly impactful genes.\n'
        '- :blue[Variance]: Ranks genes based on the variance of their loadings across components. Use this method to find genes that show '
        'varying importance across different components, capturing diverse aspects of gene expression.')
    return text


def describe_method_sb():
    text = (
        ':blue[Method of Analyzing]: Choose a :red[method] for analyzing the dataset. This helps in identifying the most relevant genes '
        'that can predict cancer outcomes such as NE, EMT, histology, grade, etc.\n'
        '\n:red[Options:]\n'
        '- :blue[Kruskal]: A non-parametric test to compare multiple groups. Use this method to identify significant differences in gene '
        'expression across various sample groups without assuming a normal distribution.\n'
        '- :blue[Neural Network]: A computational model inspired by the human brain. This method is useful for capturing complex non-linear '
        'relationships in gene expression data, making it ideal for predicting cancer outcomes based on intricate patterns.\n'
        '- :blue[Nearest Centroid]: Classifies data based on the nearest centroid. This method is effective for classifying samples into distinct '
        'groups based on gene expression profiles, providing clear predictions for different cancer types or stages.\n'
        '- :blue[Nearest Neighbor]: Classifies data based on the nearest k neighbors. Use this method to identify gene expression patterns that are '
        'most similar to those in other samples, useful for predicting outcomes by comparison to known cases.\n'
        '- :blue[Naive Bayes]: A probabilistic classifier based on Bayes\' theorem. This method assumes independence between features and is useful '
        'for rapid classification based on the probability of gene expression patterns, making it effective for initial screening and analysis.')
    return text


def describe_p_correction_kruskal():
    text = (':blue[P Value Correction Method]: Select the method for :red[correcting p-values] in multiple testing.\n'
            '\n:red[Options:]\n'
            '- :blue[fdr_bh]: False discovery rate control using Benjamini-Hochberg procedure.\n'
            '- :blue[fdr_by]: False discovery rate control using Benjamini-Yekutieli procedure.\n'
            '- :blue[bonferroni]: Bonferroni correction for multiple comparisons.\n'
            '- :blue[holm]: Holm-Bonferroni method for controlling family-wise error rate.')
    return text


def describe_ds_select2():
    text = (':blue[Choose a Database]: Select a :red[database] to serve as the basis for the analysis.\n'
            '\n:red[Options:]\n'
            '- :blue[George-SCLC]: Data from George et al. on Small Cell Lung Cancer.\n'
            '- :blue[Anish-SCLC]: Data from Anish et al. on Small Cell Lung Cancer.\n'
            '- :blue[Jiang-SCLC]: Data from Jiang et al. on Small Cell Lung Cancer.\n'
            '- :blue[George-LCNEC]: Data from George et al. on Large Cell Neuroendocrine Carcinoma.\n'
            '- :blue[Fernandez-Carcinoid]: Data from Fernandez et al. on Carcinoid Tumors.\n'
            '- :blue[Alcala-Carcinoid]: Data from Alcala et al. on Carcinoid Tumors.\n'
            '- :blue[Liu-SCLC]: Data from Liu et al. on Small Cell Lung Cancer.\n'
            '- :blue[Rousseaux-Mixed(non-NE)]: Data from Rousseaux et al. on mixed non-neuroendocrine tumors.\n'
            '- :blue[TCGA]: Data from The Cancer Genome Atlas on various tumor types.')
    return text


def describe_method_sb_ana_opt():
    text = (
        ':blue[Select a Method of Analysing]: Choose a :red[method] for analyzing the dataset. This selection helps in determining the approach to identify '
        'the most relevant genes that can predict cancer outcomes such as NE, EMT, histology, grade, etc.\n'
        '\n:red[Options:]\n'
        '- :blue[Analysing the Dataset]: Use this option to perform comprehensive analysis on the entire dataset. This method leverages various statistical '
        'and machine learning techniques to identify significant genes across all samples, providing a broad understanding of the gene expression patterns '
        'and their implications in cancer prediction.\n'
        '- :blue[Clinical Detail Based Analysing]: Select this method to focus on gene expression analysis based on specific clinical details. This approach '
        'integrates clinical data with gene expression profiles to identify genes that are significantly associated with clinical characteristics such as patient '
        'survival, disease stage, or treatment response. It helps in uncovering clinically relevant biomarkers and provides insights into the molecular mechanisms '
        'underlying different clinical phenotypes.')
    return text


def describe_clin_detail():
    text = ':blue[Clinical Detail]: Choose a :red[clinical detail] for the analysis.\n'
    return text


def describe_random_state_gene_sign2():
    text = (':blue[Random State of Your Work]: Set the :red[random state] for reproducibility in analyses. '
            'This ensures that the results are :red[consistent across multiple runs].')
    return text


def describe_num_com_ni_std():
    text = (
        ':blue[Number of Most Variable Genes]: Select the number of :red[most variable genes] based on :red[standard deviation] '
        'for further analysis. This helps in focusing on the most :red[informative genes].')
    return text


def describe_num_com_ni2():
    text = (
        ':blue[Number of Chosen Genes]: Specify the number of :red[top-ranked genes] to select for further analysis.')
    return text


def describe_n_genes_nir():
    text = ':blue[Number of Genes]: Specify the number of :red[genes] to consider in the analysis.'
    return text


def describe_var_chb2():
    text = (':blue[Use Variance Ranking]: Checkbox to enable :red[variance ranking] for gene selection. '
            'This helps in selecting genes based on their :red[variance in the dataset].')
    return text


def describe_group_ms():
    text = (':blue[Clinical Features]: Multiselect to choose :red[clinical features] for the analysis. '
            'This allows selecting multiple :red[clinical attributes] to include in the analysis.')
    return text


def describe_var_percentile2():
    text = (':blue[Variance Percentile]: Specify the :red[percentile] for :red[variance filtering]. '
            'This helps in selecting the most :red[variable genes] based on the specified threshold.')
    return text


def describe_pca_compo_slider2():
    text = (
        ':blue[Number of Components]: Select the :red[number of components] for dimensionality reduction techniques '
        'like :red[PCA, ICA, FA, SVD], etc. This determines the number of :red[principal components] to retain.')
    return text


def describe_gl_toggle2():
    text = ':blue[Use Gene List]: Toggle to use a :red[predefined gene list] for the analysis.'
    return text


def describe_n_neighbors():
    text = (
        ':blue[Number of K Neighbour]: Select the number of :red[neighbors] for the :red[K-Nearest Neighbors] classifier.')
    return text


def describe_kruskal_alpha():
    text = (':blue[Alpha Value]: Select the :red[alpha value] for the Kruskal-Wallis test. '
            'This is the :red[significance level] used for hypothesis testing.')
    return text


def cohort_description():
    text = ("\n :red[Neuroendocrine Lung Tumor Cohorts]\n"
            "\n :blue[George-SCLC:]"
            "\n81 SCLC samples, mostly treatment-naive, enriched for earlier stages. Study by George et al., published in *Nature*, 2015. [DOI: 10.1038/nature14664](https://doi.org/10.1038/nature14664)\n"
            "\n :blue[Anish-SCLC:]"
            "\n100 metastatic SCLC samples, including both NE and non-NE subtypes. Study by Lissa et al., published in *Nature Communications*, 2022. [DOI: 10.1038/s41467-022-29289-1](https://www.nature.com/articles/s41467-022-29517-9)\n"
            "\n :blue[Jiang-SCLC:]"
            "\n50 Chinese SCLC samples, primarily treatment-naive. Study by Jiang et al., published in *PLOS Genetics*, 2016. [DOI: 10.1371/journal.pgen.1005895](https://doi.org/10.1371/journal.pgen.1005895)\n"
            "\n :blue[Liu-SCLC:]"
            "\n107 SCLC samples, mostly chemo-naive. Study by Liu et al., published in *Cell*, 2023. [DOI: 10.1016/j.cell.2023.12.004](https://doi.org/10.1016/j.cell.2023.12.004)\n"
            "\n :blue[George-LCNEC:]"
            "\n66 LCNEC samples, with two molecular subtypes identified. Study by George et al., published in *Nature Communications*, 2018. [DOI: 10.1038/s41467-018-03099-x](https://doi.org/10.1038/s41467-018-03099-x)\n"
            "\n :blue[Fernandez-Carcinoid:]"
            "\n65 pulmonary carcinoid samples, with frequent mutations in chromatin-remodelling genes. Study by Fernandez-Cuesta et al., published in *Nature Communications*, 2014. [DOI: 10.1038/ncomms4518](https://doi.org/10.1038/ncomms4518)\n"
            "\n :blue[Alcala-Carcinoid:]"
            "\n115 pulmonary carcinoid samples, including 37 atypical cases, characterized by integrative genomic analysis. Study by Alcala et al., published in *Nature Communications*, 2019. [DOI: 10.1038/s41467-019-11276-9](https://doi.org/10.1038/s41467-019-11276-9)\n"
            "\n :blue[Rousseaux-Mixed Lung Tumor:]"
            "\nVarious types of lung cancer, including ADC, SQC, BAS, SCLC, LCNEC, carcinoid, LCC, and unidentified tumors, totaling 307 samples. Study by Rousseaux et al., published in *Nature Genetics*, 2013. [DOI: 10.1126/scitranslmed.3005723](https://www.science.org/doi/10.1126/scitranslmed.3005723)\n"
            "\n :red[* For more information of subtypes of cohorts, please visit 'Subtype Selector' TAB *]\n"
            "\n :blue[* Data Preprocessing Information *]"
            "\n :red[* Only samples with transcriptomic data were included in the analysis. For survival analysis, only samples with both transcriptomic and survival data were used. *]\n"
            "\n :red[* Gene names were standardized to HUGO symbols if the original study used a different nomenclature. Genes that could not be precisely mapped to HUGO symbols were removed from the dataset. *]\n"
            "\n :red[* Quality control was applied to ensure that samples with low RNA integrity or other quality issues were excluded from analysis. *]\n"
            "\n :red[* Batch effects were assessed and corrected using appropriate normalization techniques to ensure comparability across datasets. *]\n"
            "\n :red[* Data was log-transformed and normalized to account for differences in sequencing depth and to stabilize variance. *]\n")

    return text

