-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Pipeline-PGPR_AgroBioTech

-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#### Objetivo

O objetivo deste pipeline é utilizar genomas anotados e depositados no Genbank, hospedado no NCBI, de bactérias que estão relacionadas à promoção de crescimento de plantas, para auxiliar na análise \textit{in silico} do genoma das cepas \textit{Bacillus velezensis} LABIMs, que são anotações genômicas disponibilizadas pela Universidade Estadual de Londrina (UEL) para que, a partir das análises realizadas, construir um banco de dados de genes relacionados a promoção de crescimento em plantas, a fim de corroborar com estudos futuros nesta área de pesquisa. 

-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#### A sequência de scripts a serem compilados neste pipeline é: 

1. download_e_organizacao_genomas.py
2. blastp.py
3. analisar_blastp.py
4. comparacao_genomas.py
5. construcao_banco_AgroBioTech.py

Observações:

- Para compilar o script 1, recomenda-se utilizar um cluster. Nesta pesquisa, foi utilizado os equipamentos do Laboratório Multiusuários CCCT-CP, que utiliza o software SLURM para o gerenciamento de filas de trabalho. Utilizamos o script shell "script_python_quali.sh" para compilar o script 1.

-------------------------------------------------------------------------------------------------------------------------------------------------------------------------


#### Bibliotecas Python e softwares necessários:

1. Python - https://www.python.org/downloads/
2. Biblioteca Biopython - https://biopython.org/wiki/Download
3. BLAST+ - https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/
   

-------------------------------------------------------------------------------------------------------------------------------------------------------------------------


### Fluxograma do Pipeline 

<div align="center">
<img src="https://github.com/Biahsilva/Pipeline_PGPR_AgroBioTech/assets/102994978/f407b473-5392-401a-a904-029381c2e6fe.png" width="600px" /> 
</div>
