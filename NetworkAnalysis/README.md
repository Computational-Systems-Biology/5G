# Network Analysis

To investigate the network coherence of the differentially expressed/methylated genes in the context of biological networks, we employed two gene-centric metabolic networks (Recon2 \cite{Thiele2013} and Recon3D \cite{Brunk2018}) and two gene-level protein-protein interaction networks (String \cite{Szklarczyk2021} and Biogrid \cite{Oughtred2021}). Our focus was on the effective subnetworks, which represent the projections of differentially expressed/methylated genes onto the employed gene-centric metabolic network or gene-level protein-protein interaction network. We conducted a comprehensive analysis of network coherence on these subnetworks, by thoroughly assessing whether the network coherence happens to be greater or lower than expected at random.

## Gene-centric metabolic network construction:
We extracted gene-centric metabolic models (GCMNs) from the Recon2 \cite{Thiele2013} and Recon3D \cite{Brunk2018} human metabolic network models following the methods described in \cite{palsson2006systems,sonnenschein2011analog,sonnenschein2012network,knecht2016distinct}. In GCMNs, nodes represent genes and edges represent the associations of genes via metabolic reactions. Stated differently, a connection between two genes is established if the metabolic reactions associated with these genes share a common metabolite.

The primary exchange metabolites, namely ATP, ADP, CO2, H, NAD, NADH, among others, are the most highly connected metabolic species that are unlikely to establish links between genes with similar metabolic functions, leading to an artificially denser metabolic network \cite{ma2003connectivity,ma2003reconstruction,kharchenko2005expression}. In order to mitigate this effect, we opt to remove a subset of metabolites that represent the top 2\% of the most highly connected metabolites, prior to network construction.

The resulting gene-centric metabolic network obtained through the utilization of the Recon2 human metabolic network model has 1806 nodes and 31699 edges. Similarly, the one extracted from Recon3D has 3449 nodes and 233235 edges.

## Gene-level protein-protein interaction network construction:
We perform our analysis using two different protein-protein interaction networks on gene-level (GPINs) derived from String \cite{Szklarczyk2021} and Biogrid \cite{Oughtred2021} protein-protein interaction databases.

To construct GPIN derived from the String database, we selectively included protein interactions associated with the human organism, and only the ones with a score exceeding 850 to include solely direct interactions. We then cross-referenced the protein IDs with associated gene IDs using the Ensembl database \cite{Cunningham2022}. The GPIN is depicted as a graph, with genes serving as nodes and the edges representing the associations between the genes through protein interactions. Notably, the graph has 12468 nodes and 139560 edges.

The GPIN derived from the Biogrid database is constructed following the methodology explained above without any filtration of interaction scores. The graph consists of a total of 21797 nodes and 1018531 edges.

## ID Mapping:
We employed the Python mygene package to map gene names to both Entrez IDs and Ensembl IDs, for the purpose of performing GO and network analyses, respectively. In instances where multiple hits are identified, meaning that a singular gene name may correspond to multiple IDs, all of the mapped IDs were included.

## Network Coherence:
The differentially expressed/methylated genes were projected onto the network employed to extract the effective subnetwork. The coherence value of a subnetwork was determined by the ratio of non-isolated genes to the total number of genes that exist within the network. The analysis was only performed if the number of differentially expressed/methylated genes was greater than or equal to 5.

To obtain the null distribution of coherence values, 5000 gene sets, each being equivalent in size to the effective subnetwork, are randomly drawn from the employed network, and the coherence values of these randomly drawn subnetworks were subsequently calculated. The z-score of the coherence of the effective subnetwork was then determined utilizing this distribution.
