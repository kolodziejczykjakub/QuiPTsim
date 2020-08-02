library(DiagrammeR)

my_graphviz <- grViz("digraph{
graph[rankdir = TD]
node[shape = rectangle, style=filled, margin = 0.35]
A[label = 'Create n-grams \\n (window size = 10)']
B[label = 'Filter out features \\n non informative k-mers']
H[label = 'H0 is rejected \\n if IG > threshold']
I[label = 'Benjamini-Hochberg \\nProcedure']
J[label = 'Filtered n-grams']
subgraph cluster_1 {
  graph[shape = rectangle]
  style = rounded
  bgcolor = Gold

  label = 'Fisher`s exact test'
  node[shape = rectangle, fillcolor = LemonChiffon, margin = 0.25]
  C[label = 'Calculating p-values']
  E[label = 'Exact p-value \\n if n11 <= 200']
  F[label = 'Monte Carlo approach \\n if n11 > 200']
  G[label = 'Information Gain']
}
edge[color = black]
A -> B
B -> C
C -> E
C -> F
E -> I
F -> I
I -> H
H -> J
}")

my_graphviz
