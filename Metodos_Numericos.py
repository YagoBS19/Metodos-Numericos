import matplotlib.pyplot as plt #Import das configurações utilizadas no programa
import os
import sys
from sympy import *
import sympy as sp
valores_x, valores_y, valores_k, valores_e, valores_xe  = [], [], [], [], [] #Definição inicial de valores
x0,y0,nav,i,j,k = 0,0,0,0,0,0
equacao_inicial = ''
x = sp.Symbol('x')#Definição de X e Y para a resolução da EDO
y = sp.Function('y')(x)

#Código de resolução de qualquer EDOHomogênea de Primeira Ordem
#O programa só necessita dos valores do lado direito da equação
#Exemplo: Ao invés de digitar "y'=y*cos(x)" ou "dy/dx=y*cos(x)" o usuário deverá digitar somente "y*cos(x)"
#Glossário: Tangente = Tan(x), Seno = Sin(x), Cosseno = Cos(x), Potência = y**x, Log = Log(x), Ln = Ln(x)
def equacao(x, y):#Definição da EDO
    equacao_final = eval(equacao_inicial)
    return equacao_final

def menu():#Função de Impressão das opções do menu
   print('0-Solução Exata')
   print('1-Método de Euler')
   print('2-Método de Euler aprimorado')
   print('3-Runge-Kutta de ordem 3')
   print('4-Runge-Kutta de ordem 4')
   print('5-Finalizar o programa\n')

def imprimir_exata(valores_e,valores_xe,j):#Função de Impressão da Solução Exata
    print('----------------------------------------------------------------------------')
    print('O valor de Y para '+str(j+1)+'º iteração é:',valores_e[j])
    print('O valor de X para a '+str(j+1)+'º iteração é: ',valores_xe[j])
    print('----------------------------------------------------------------------------')

def imprimir_simples(valores_x,valores_y,j):#Função de impressão de Runge-Kutta Ordem 1 e 2
    print('----------------------------------------------------------------------------')
    print('O valor de Y para '+str(j+1)+'º iteração é:',valores_y[j])
    print('O valor de X para a '+str(j+1)+'º iteração é: ',valores_x[j])
    print('----------------------------------------------------------------------------')

def imprimir(valores_x,valores_y,valores_k,j):#Função de impressão de Runge-Kutta Ordem 3 e 4
    print('----------------------------------------------------------------------------')
    print('O valor de Y para '+str(j+1)+'º iteração é:',valores_y[j])
    print('O valor de X para a '+str(j+1)+'º iteração é: ',valores_x[j])
    print('Os valores de K para a '+str(j+1)+'º iteração é:',valores_k)
    print('----------------------------------------------------------------------------')

def imprimir_erro(valores_y,valores_e,j):#Função de impressão de erro
    print('----------------------------------------------------------------------------')
    print('O valor de Y para '+str(j+1)+'º iteração é:',valores_y[j])
    print('O valor exato para a '+str(j+1)+'º iteração é: ',valores_e[j])
    print('O erro correspondente a '+str(j+1)+'º iteração é:',(valores_e[j]-valores_y[j]))
    print('----------------------------------------------------------------------------')

def Exata(x0,y0,x,y,i,j):#Função de resolução da Solução Exata
    equacao_diferencial = sp.Eq(y.diff(x), equacao(x,y))
    ics = {y.subs(x, x0): y0}
    solucao_geral = sp.dsolve(equacao_diferencial, ics=ics)
    x0 = round(x0+h,2)
    while j!=i:
        equacao_exata = solucao_geral.subs(x, x0)
        valores_e.insert(j, equacao_exata.rhs)
        valores_xe.insert(j, x0)
        x0 = round(x0+h,2)
        j = j+1

def Ordem1(x0,y0,i,j):#Função de resolução da Ordem 1
    while j!=i:
        yi = y0+equacao(x0,y0)*h
        valores_y.insert(j,yi)
        valores_x.insert(j,x0)
        imprimir_simples(valores_x,valores_y,j)
        y0 = yi
        x0=round(x0+h,2)
        j=j+1

def Ordem2(x0,y0,i,j):#Função de resolução da Ordem 2
    while j!=i:
        Euler = y0+equacao(x0,y0)*h
        yi = y0+((equacao(x0,y0)+equacao(x0+h,Euler))/2)*h
        valores_y.insert(j,yi)
        valores_x.insert(j,x0)
        imprimir_simples(valores_x,valores_y,j)
        y0 = yi
        x0=round(x0+h,2)
        j=j+1
          
def Ordem3(x0,y0,i,j):#Função de resolução da Ordem 3
    while j!=i:
        k1 = equacao(x0,y0)*h
        k2 = equacao(x0+h/2,y0+k1/2)*h
        k3 = equacao(x0+3*h/4,y0+3*k2/4)*h
        yi = y0+((2*k1+3*k2+4*k3)/9)
        valores_y.insert(j,yi)
        valores_x.insert(j,x0)
        valores_k = [k1,k2,k3]
        imprimir(valores_x,valores_y,valores_k,j)
        valores_k = []
        y0 = yi
        x0 = round(x0+h,2)
        j = j+1
    
def Ordem4(x0,y0,i,j):#Função de resolução da Ordem 4
    while j!=i: 
        k1 = equacao(x0,y0)*h
        k2 = equacao(x0+h/2,y0+k1/2)*h
        k3 = equacao(x0+h/2,y0+k2/2)*h
        k4 = equacao(x0+h,y0+k3)*h
        yi = y0+((k1+2*k2+2*k3+k4)/6)
        valores_y.insert(j,yi)
        valores_x.insert(j,x0)
        valores_k = [k1,k2,k3,k4]
        imprimir(valores_x,valores_y,valores_k,j)
        valores_k = []
        y0=yi
        x0=round(x0+h,2)
        j = j+1

def navigation(nav,y,x0,y0): #Função do menu de seleção de resolução
    match nav:
        case 0:
            j=0
            Exata(x0,y0,x,y,i,j)
            while j!=i:
                imprimir_exata(valores_e,valores_xe,j)
                j=j+1
        case 1:
            j=0
            print('-----Tabela de valores comparativos-----\n') 
            Exata(x0,y0,x,y,i,j)
            Ordem1(x0,y0,i,j)
            os.system('pause')
            os.system('cls')
            print('-----Tabela de erros-----\n')
            while j!=i:
                imprimir_erro(valores_y,valores_e,j)
                j=j+1
        case 2:
            j=0
            Exata(x0,y0,x,y,i,j)
            Ordem2(x0,y0,i,j)
            os.system('pause')
            os.system('cls')
            print('-----Tabela de erros-----\n')
            while j!=i:
                imprimir_erro(valores_y,valores_e,j)
                j=j+1
        case 3:
            j=0
            print('-----Tabela de valores comparativos-----\n')
            Exata(x0,y0,x,y,i,j)
            Ordem3(x0,y0,i,j)
            os.system('pause')
            os.system('cls')
            print('-----Tabela de erros-----\n')
            while j!=i:
                imprimir_erro(valores_y,valores_e,j)
                j=j+1
        case 4:
            j=0
            print('-----Tabela de valores comparativos-----\n')
            Exata(x0,y0,x,y,i,j)
            Ordem4(x0,y0,i,j)
            os.system('pause')
            os.system('cls')
            print('-----Tabela de erros-----\n')
            while j!=i:
                imprimir_erro(valores_y,valores_e,j)
                j=j+1

while nav!=5:#Loop para o menu

    #Início da configuração de loop para o menu
    os.system('cls')
    menu()
    nav = int(input('Informe a ordem que deseja utilizar: '))
    if nav == 5:
        os.system('cls')
        print('Obrigado por utilizar o programa!')
        sys.exit(0)
    while nav < 0 or nav > 5:
       os.system('cls')
       print('Ordem inválida.')
       os.system('pause')
       os.system('cls')
       menu()
       nav = int(input('Informe a ordem que deseja utilizar: '))
       if nav == 5:
          print('Obrigado por utilizar o programa!')
          sys.exit(0)
    os.system('cls')
    equacao_inicial = str(input('Informe uma EDO homogênea de primeiro grau: '))
    x0 = float(input('Informe o valor de X0: '))
    y0 = float(input('Informe o valor de Y0: '))
    h = float(input('Informe o valor de h: '))
    i = int(input('Informe o número de iterações: '))
    os.system('cls')
    #Fim da configuração de loop para o menu

    valores_x, valores_y, valores_e, valores_k, valores_xe = [], [], [], [], []

    navigation(nav,y,x0,y0) #Menu para a seleção de resolução

    if valores_y != []: #Plotagem de gráfico para os Métodos Numéricos
        fig, ax = plt.subplots()
        ax.plot(valores_x, valores_y, color='blue')
        ax.plot(valores_xe, valores_e, color='red')
        ax.set(xlabel='X', ylabel='Y')
        ax.grid()
        plt.suptitle('Gráfico de Soluções')
        plt.title('Azul: Métodos Numéricos e Laranja: Solução Exata')
        plt.show()
        ax.cla()

    else: #Plotagem de gráfico para a Solução Exata
        fig, ax = plt.subplots()
        ax.plot(valores_xe, valores_e, color='red')
        ax.set(xlabel='X', ylabel='Y')
        ax.grid()
        plt.suptitle('Gráfico de Soluções')
        plt.show()
        ax.cla()
    os.system('pause')