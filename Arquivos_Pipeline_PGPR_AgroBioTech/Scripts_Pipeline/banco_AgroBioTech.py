from Bio import SeqIO, Entrez
from Bio.Blast.Applications import*
from Bio.Blast import NCBIXML
import os
import sys
import csv

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

def criar_banco_AgroBioTech():

    arquivo_referencia = open('../banco_de_dados_AgroBioTech/genes_pgp_referencias/genes_pgp_referencias_AgroBioTech.fasta','r')
    arquivo_bs96 = open('../banco_de_dados_AgroBioTech/genes_pgp_adicionais/genes_pgp_CP090125.1.fasta','r')
    arquivo_rz2ms9 = open('../banco_de_dados_AgroBioTech/genes_pgp_adicionais/genes_pgp_CP049978.1.fasta','r')

    with open("../Banco_de_dados/banco_de_dados_AgroBioTech/banco_de_genes_AgroBioTech/genes_pgp_AgroBioTech.fasta","w") as arquivo_AgroBioTech:
        
        arquivo_AgroBioTech.write(arquivo_referencia.read())
        arquivo_AgroBioTech.write(arquivo_bs96.read())
        arquivo_AgroBioTech.write(arquivo_rz2ms9.read())

    arquivo_AgroBioTech.close()
    arquivo_referencia.close()
    arquivo_bs96.close()
    arquivo_rz2ms9.close()

def download_queries(genoma_id):

    saida = open('../Banco_de_dados/banco_de_dados_adicionais/dados_genes_sbjct_'+str(genoma_id)+'.fasta','w')
    handle = Entrez.efetch(db='nuccore', id=genoma_id, rettype='gb')

    for seq_record in SeqIO.parse(handle, "genbank") :
        
        for feature in seq_record.features:
            
            if feature.type=="CDS" and 'translation' in feature.qualifiers.keys():

                quali = ""

                for keys in feature.qualifiers.keys():

                    if not keys == 'translation':
                        quali += " ".join(feature.qualifiers[keys])+' '

                if 'gene' in feature.qualifiers.keys():

                    saida.write(">%s %s\n%s\n" % (
                    feature.qualifiers['gene'][0],
                    quali,
                    feature.qualifiers['translation'][0]))

                else:

                    saida.write(">%s\n%s\n" % (
                    quali,
                    feature.qualifiers['translation'][0]))

    saida.close()
    handle.close()

def blastp(genoma_id, similaridade):

    query = "../Banco_de_dados/banco_de_dados_AgroBioTech/banco_de_genes_AgroBioTech/genes_pgp_AgroBioTech.fasta"
    subject = '../Banco_de_dados/banco_de_dados_adicionais/dados_genes_sbjct_'+str(genoma_id)+'.fasta'
    comando_blastp = NcbiblastpCommandline(query=query, subject=subject, max_target_seqs=1, out="../Relatorios_blast/relatorios_blast_comparacao/"+similaridade+"/relatorio_blastp_"+str(genoma_id)+".xml", outfmt=5)
    stdout, stderr = comando_blastp()

def ler_blastp(genoma_id, similaridade):

    aux = 1
    result_handle = open("../Relatorios_blast/relatorios_blast_comparacao/"+similaridade+"/relatorio_blastp_"+str(genoma_id)+".xml", 'r')
    blast_iterator = NCBIXML.parse(result_handle)

    if os.path.exists("../Relatorios_blast/relatorios_blast_final/"+similaridade+"/novo_relatorio_blastp_"+str(genoma_id)+".txt") and aux==1:
            os.remove("../Relatorios_blast/relatorios_blast_final/"+similaridade+"/novo_relatorio_blastp_"+str(genoma_id)+".txt")
    
    arquivo_final = open("../Relatorios_blast/relatorios_blast_final/"+similaridade+"/novo_relatorio_blastp_"+str(genoma_id)+".txt",'a')
    arquivo_adicional = open('../Resultado_similaridade/'+similaridade+'/resultado_similaridade_'+str(genoma_id)+'.csv', 'a')
    writer = csv.writer(arquivo_adicional)
    header = ['gene', 'presenca']
    writer.writerow(header)
    i = 0

    for blast_record in blast_iterator:

        if not blast_record.alignments:
            
            data = [str(blast_record.query).split(' ', 1)[0], '1']
            writer.writerow(data)

        else:

            hsp = blast_record.alignments[0].hsps[0]
            alignment = blast_record.alignments[0]

            if hsp.identities/alignment.length > 0.9:

                aux = 0


                arquivo_final.write('ALINHAMENTO BLAST '+str(blast_record.query)+'\n\n')

                arquivo_final.write('Sequencia: '+str(alignment.title)+'\n')
                arquivo_final.write('Tamanho: '+str(alignment.length)+'\n')
                arquivo_final.write('Score: '+str(hsp.score)+'\n')
                arquivo_final.write('Bits: '+str(hsp.bits)+'\n')
                percentage = "{:.0%}".format(hsp.identities/alignment.length)
                arquivo_final.write('Identities: '+str(percentage)+'\n')
                arquivo_final.write('E-Value: '+str(hsp.expect)+'\n')
                arquivo_final.write(str(hsp.query)+'\n')
                arquivo_final.write(str(hsp.match)+'\n')
                arquivo_final.write(str(hsp.sbjct)+'\n')
                arquivo_final.write("\n\n\n")

                data = [str(blast_record.query).split(' ', 1)[0], '0']

            else:

                data = [str(blast_record.query).split(' ', 1)[0], '1']

            writer.writerow(data)

    
    arquivo_final.close()
    arquivo_adicional.close()

def analisar(lista_id, similaridade):

    dicionario = {}

    for genoma_id in lista_id:

        with open('../Resultado_similaridade/'+similaridade+'/resultado_similaridade_'+str(genoma_id)+'.csv') as arquivo_csv:
            
            leitor_csv = csv.reader(arquivo_csv, delimiter=',')

            for linhas in leitor_csv:

                if linhas[1] == '0':

                    if linhas[0] in dicionario.keys():

                        dicionario[linhas[0]] += 1

                    else:

                        dicionario[linhas[0]] = 1

                elif linhas[1] == '1':

                    if linhas[0] not in dicionario.keys():

                        dicionario[linhas[0]] = 0
    
    arquivo = open('../Resultado_similaridade/'+relatorio+'/analise_similaridade.txt','w')
    arquivo.write('\n'.join(str(valor) for valor in [chave for chave, valor in dicionario.items() if valor == 11]))
    arquivo.close()

# Leitura do Input do Script
if len(sys.argv) != 3:

    print("Uso: python blast_similaridade.py input.txt similaridade")

else:

    nome_arquivo = sys.argv[1]
    similaridade = str(sys.argv[2])
    valores = ler_arquivo(nome_arquivo)

    if valores:

        criar_banco_AgroBioTech()

        for genoma_id in valores:

            download_queries(genoma_id)
            blastp(genoma_id, similaridade)
            ler_blastp(genoma_id, similaridade)

        analisar(valores, similaridade)

#Fim do Script
