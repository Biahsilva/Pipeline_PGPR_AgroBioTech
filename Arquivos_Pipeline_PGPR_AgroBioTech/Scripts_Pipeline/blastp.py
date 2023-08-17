from Bio import Entrez, SeqIO
from Bio.Blast.Applications import*
from Bio.Blast import NCBIXML
import os
import sys

Entrez.email = 'beatrizsilva@alunos.utfpr.edu.br'
Entrez.api_key = '3d817b61bb5164b63368b92a9f158f8d4308'

# Função para ler o arquivo Input
def ler_arquivo(nome_arquivo):
    try:
        with open(nome_arquivo, 'r') as arquivo:
            linhas = arquivo.readlines()
            valores = [linha.strip() for linha in linhas]
            return valores
    except FileNotFoundError:
        print("Arquivo não encontrado.")
        return []

def download_queries(genoma_id):

    arquivo = open('../Codigos_busca_genes/Listas_genes_pgp_faltantes/codigos_genes_pgp_faltantes_'+str(genoma_id)+'.txt','r')
    aux = 1

    genes = arquivo.readlines()

    for gene in genes:

        gene = gene.strip()

        if os.path.exists('../Banco_de_dados/banco_de_dados_blast/banco_de_dados_queries/query_'+str(genoma_id)+'.fasta') and aux == 1:
            os.remove('../Banco_de_dados/banco_de_dados_blast/banco_de_dados_queries/query_'+str(genoma_id)+'.fasta')

        saida = open('../Banco_de_dados/banco_de_dados_blast/banco_de_dados_queries/query_'+str(genoma_id)+'.fasta','a')
        handle = Entrez.efetch(db='protein', id=gene, rettype='fasta')
        seqRecord = SeqIO.read(handle, format='fasta')
        handle.close()
        saida.write('>'+str(seqRecord.description)+'\n')
        saida.write(str(seqRecord.seq)+'\n')
        saida.close()
        aux = 0

def blastp(genoma_id):

    query = "../Banco_de_dados/banco_de_dados_blast/banco_de_dados_queries/query_"+str(genoma_id)+".fasta"
    subject = "../Banco_de_dados/banco_de_dados_blast/banco_de_dados_sbjcts/dados_genes_sbjct_"+str(genoma_id)+".fasta"
    comando_blastp = NcbiblastpCommandline(query=query, subject=subject, max_target_seqs=1, out="../Relatorios_blast/relatorios_blast_complementar/relatorio_blastp_"+str(genoma_id)+".xml", outfmt=5)
    stdout, stderr = comando_blastp()

# Leitura do Input do Script
if len(sys.argv) != 2:

    print("Uso: python blastp.py input.txt")

else:

    nome_arquivo = sys.argv[1]
    valores = ler_arquivo(nome_arquivo)

    if valores:

        for genoma_id in valores:

            download_queries(genoma_id)
            blastp(genoma_id)

#FIM
