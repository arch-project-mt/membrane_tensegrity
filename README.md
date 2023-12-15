# membrane_tensegrity
## Initial setup
- Install Docker ([Official page](https://docs.docker.com/compose/install/))

- You can check if docker and docker-compose are available on your computer by the following command

```bash
docker -v
```
- If the message like the following is returned, it means you successfully installed Docker
```bash
Docker version 20.10.13, build a224086
```
- Check if docker-compose is available on your computer in the same way
```bash
docker-compose -v
```

- After the docker and docker-compose are set up, you can execute the following command to activate the docker environment on your computer

```bash
docker-compose up -d
```
- You can enter the docker container by the following command
```bash
docker container exec it mt bash
```
- After entering the docker container, you can compile the C++ file for computing RMSD

```bash
bash build_rmsd.sh
```

## Project workflow
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
## TODO
- Develop and implement the algorithm for solving the following optimization problem

$$
\begin{align}
&\min_{M_s} RMSD(M_t, M_s),\notag\\
&\text{where }RMSD(M_t, M_s) \notag\\
&= \min_{\mathbf{R}, \vec{v}}\sqrt{\frac{1}{n}\sum_{i=1}^n||\vec{m_{t,i}} - \mathbf{R}(\vec{m_{s,i}} - \vec{v})||^2},\notag\\
&s.t. \text{ developability constraints}\notag
\end{align}
$$
