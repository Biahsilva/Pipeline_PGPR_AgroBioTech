import csv
import sys

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

def analisar(lista_id):

    relatorio = '90'
    dicionario = {}

    for genoma_id in lista_id:

        with open('../Resultado_similaridade/'+relatorio+'/resultado_similaridade_'+str(genoma_id)+'.csv') as arquivo_csv:
            
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
if len(sys.argv) != 2:

    print("Uso: python download_e_organizacao_genomas.py input.txt")

else:

    nome_arquivo = sys.argv[1]
    valores = ler_arquivo(nome_arquivo)

    if valores:

        criar_pastas()

        for genoma_id in valores:

            download_genoma(genoma_id)
            organizar_dados(genes_genoma_id)
            buscar_genes(genes_genoma_id)

#Fim do Script


analisar()

#FIM
