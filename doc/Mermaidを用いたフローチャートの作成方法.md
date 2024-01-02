# Mermaidを用いたフローチャートの作成方法

Created by: KOYANOBunsho
Created time: December 16, 2023 1:03 PM
Tags: Markdown

## 概要

- Groupworkの円滑な進行や進捗報告を行うためにMermaidを用いたフローチャートが便利になると考えられる
- 本ドキュメントでは，Groupworkで必要になるフローチャートを作成するためのMermaidコードをまとめる
- Mermaidコードの作成は，以下のWebページを参考にした

[Mermaid入門](https://zenn.dev/kento_mm_ninw/articles/8b10afdbef306a)

## Workflow for the whole project

```mermaid
flowchart LR
subgraph main_project
A[Decide the shape of \n membrane tensegrity] --> B[By grasshopper script,\n make 2D strut pattern] --> C[Construct the target membrane tensegrity \n with the real world materials]
end

subgraph side_project
B -->|In parallel| E[Improve the optimization process]
E --> F[Propose the new optimization problem]
E --> G[Evaluate the optimized result]
end
```

## Computational flow for the subproject

```mermaid
flowchart TD

subgraph Inside the grasshopper
B[Compute the target meshes\n and extract information on vertices and edges] -->|Input coordinates\n and adjacency lists representing edges| A[Python script]
A --> |Input the vertices and \nedges of simulated meshes| E[Complete the 2D strut patter of the simulated meshes]
end

subgraph Outside the grasshopper
A -->|Output coordinates\n and adjacency lists representing edges| C[Store the coordinates and adjacensy lists as csv files]
C --> D[Solve the optimization problem by C++]
D -->|Input the computation result| A
end
```

## Schedule

```mermaid
flowchart LR
subgraph Schedule
A --> B
end

```
