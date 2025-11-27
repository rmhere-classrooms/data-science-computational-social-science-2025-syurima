# II Grafy preferential attachment (Barabási-Albert)

library(igraph)

# Barabasi-Albert graph with 1000 nodes
g_ba <- barabasi.game(n=1000)

# Visualize graph with Fruchterman & Reingold layout
layout <- layout.fruchterman.reingold(g_ba)
plot(g_ba, layout=layout, vertex.size=2,
     vertex.label=NA, edge.arrow.size=.2)

# Find the most central node according to betweenness
betweenness_values <- betweenness(g_ba)
most_central_node <- which.max(betweenness_values)
most_central_node
# [1] 13 -> najbardziej centralny węzeł to węzeł nr. 13

# Calculate the diameter of the graph
graph_diameter <- diameter(g_ba)
graph_diameter
# [1] 9 -> średnica grafu wynosi 9


# ---------------------------------
# 6. W komentarzu napisz czym różnią się grafy Barabási-Albert i Erdős-Rényi.

# Grafy te różnią się sposobem generowania krawędzi między węzłami. 

# W modelu Erdős-Rényi krawędzie są tworzone losowo z jednakowym prawdopodobieństwem dla każdej pary węzłów
# (co prowadzi do równomiernego rozkładu stopni węzłów - nie ma wyraźnych hubów).

# W modelu Barabási-Albert krawędzie są tworzone na zasadzie "preferencyjnego przyłączania" - nowe węzły 
# mają większe prawdopodobieństwo połączenia się z już istniejącymi węzłami o wyższym stopniu. 
# Powstaje w ten sposób sieć, gdzie kilka węzłów ma bardzo wysoki stopień (huby), ale większość ma stopień dosyć niski
# Wizualnie - rozgałęziona sieć z kilkoma wyraźnymi centrami. Węzły na gałęziach łączą się z hubami ale nie ze sobą nawzajem.
# ---------------------------------