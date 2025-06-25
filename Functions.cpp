#include <iostream>
#include <iomanip>
#include <fstream>
#include <algorithm>
#include <math.h>  
#include <vector>
#include <list>
#include <iterator>   
#include <stdlib.h>     
#include <time.h>       
#include <chrono>   
#include "Functions.h"
#include <stdarg.h>
#include <values.h>
#include <string.h>
#include <malloc.h>
#include <ctime>
#include <ilcplex/ilocplex.h>

using namespace std;

void Read_instances(char *argv, int &m, int &n, vector<double> &p, vector<double> &f, vector<double> &d,  vector<vector<double>> &c,  vector<double> &pen, double &demanda_total, double &capacidade_total, bool &coeficientes_inteiros) 
    {    
        int AUX, id, AUX2;
        ifstream arq(argv);
		if (!arq) {cerr << "Erro arquivo \n"; exit(0);}
        cout << "Instance: " << argv << endl;
        arq >> AUX;
        arq >> AUX2;
        if (AUX == 100) //beasley instances (100 locations)
        {
            id = 1;
            m = AUX;
            n = AUX2;
        }
        if (AUX > 100) //TBED 1 instances  ( > 100 locations)
        {
            id = 2;
            n = AUX;
            m = AUX2;
        }
        if ( (AUX < 100) && (AUX2 > AUX) ) //yang instances //gaadegard //holmberg  ( < 100 locations)
        {
            id = 3;
            m = AUX;
            n = AUX2;
        }
        if ( (AUX < 100) && (AUX2 < AUX) ) //Diaz instances
        {
            id = 4;
            n = AUX;
            m = AUX2;
        }
            p.resize(m);
            f.resize(m);
            d.resize(n);
            c.resize(n, vector<double>(m));
            pen.resize(m); // Caso seja usado

       if (id == 1)  //or 4
       {

        vector<vector<double>> caux(n, vector<double>(m)); // declara um vetor de vetores para representar uma matriz com n linhas e m colunas para ler os custos total de n para m
        for (int j = 0; j <  m; j++){
        arq >> p[j]; //lê as capacidades
        capacidade_total = capacidade_total + p[j];
        arq >> f[j]; //lê os custos fixos
        }
for (int i = 0; i <  n; i++){
arq >> d[i];      //lê a demanda do cliente n
demanda_total = demanda_total + d[i];
        for (int j = 0; j < m; j ++)
        {
        arq >> caux[i][j] ;
        }
        }
        for (int i = 0; i <  n; i++){
                for (int j = 0; j <  m; j++){

                c[i][j] = (double) caux[i][j]/d[i];
                }


        }
 
       }
       if (id == 2)   {
      vector<vector<double>> caux(m, vector<double>(n)); // declara um vetor de vetores para representar uma matriz com n linhas e m colunas para ler os custos total de n para m
    for (int i = 0; i <  n; i++){
        arq >> d[i];
        demanda_total = demanda_total + d[i];
        }
    for (int j = 0; j <  m; j++){
arq >> p[j];
     capacidade_total = capacidade_total + p[j];
        }
     for (int j = 0; j <  m; j++){
arq >> f[j];
        }
     for (int j = 0; j < m; j ++){
   for (int i = 0; i <  n; i++)
        {
arq >> caux[j][i];
        }
    }
     for (int j = 0; j < m; j ++){  for (int i = 0; i <  n; i++)  {  c[i][j] = caux[j][i]; } }
   
    }
        if (id == 3)
   {
        coeficientes_inteiros = true;      

            for (int j = 0; j <  m; j++)
            {
                arq >> p[j]; //lê as capacidades
                     capacidade_total = capacidade_total + p[j];
                arq >> f[j]; //lê os custos fixos
            }
            for (int i = 0; i <  n; i++){
                arq >> d[i];
                demanda_total = demanda_total + d[i];
            }
       
        vector<vector<double>> caux(m, vector<double>(n)); // 
        for (int j = 0; j < m; j ++){
            for (int i = 0; i <  n; i++)
            {
arq >> caux[j][i];
        }
    }
        for (int j = 0; j < m; j ++){
   for (int i = 0; i <  n; i++)
        {

             c[i][j] = (double) caux[j][i]/d[i];
        }
    }
    }
if (id == 4)
{
               coeficientes_inteiros = true;  
   for (int i = 0; i <  n; i++){
       for (int j = 0; j < m; j ++)
        {
            arq >> c[i][j];
        }   
    }
    for (int i = 0; i < n; i++)
    {
       arq >> d[i];
    }
        for (int j = 0; j < m; j++)
    {
       arq >> f[j];
    }
            for (int j = 0; j < m; j++)
    {
       arq >> p[j];
    }
 }
    }

double tempo_decorrido(clock_t& inicio_CPU)
{
    clock_t fim_CPU = clock();
    return (double)(fim_CPU - inicio_CPU)/CLOCKS_PER_SEC;
}

void CIgreedyA(vector<int> vetor_N1, vector<int> vetor_N2, vector<double> aux_denso, double &soma_demanda_cobertura, double B, vector<double> d, vector<int> &vetor_cobertura, bool &corte_encontrado, vector<double> x, double y_til,  vector<int> &best_CI, double &valor_best_CI)
{
double soma_demandas_adicionadas = 0;  
funcao_ordena_decrescente2(vetor_N1,aux_denso, d);
for (int i = 0; i < vetor_N1.size(); i++)
{
   if (soma_demanda_cobertura <=  B  )
   {
       vetor_cobertura.push_back(vetor_N1[i]);
       soma_demandas_adicionadas = soma_demandas_adicionadas + d[vetor_N1[i]];
       soma_demanda_cobertura = soma_demanda_cobertura + d[vetor_N1[i]];
   }
   else
   {

       break;
   }
}
if (soma_demanda_cobertura <= B)
{
funcao_ordena_decrescente(vetor_N2, d);
for (int i = 0; i < vetor_N2.size(); i++)
{
   if (soma_demanda_cobertura <=  B )
   {
       vetor_cobertura.push_back(vetor_N2[i]);
       soma_demandas_adicionadas = soma_demandas_adicionadas + d[vetor_N2[i]];
       soma_demanda_cobertura = soma_demanda_cobertura + d[vetor_N2[i]];
   }
   else
   {
       break;
   }
}
}


double soma_aux1 = 0;
for (int i = 0; i < vetor_cobertura.size(); i++)
{
     soma_aux1 = soma_aux1 + x[vetor_cobertura[i]];
}
double limitante_minimal = vetor_cobertura.size() - 1 ;
if (soma_aux1 - 0.000000001 > limitante_minimal*y_til)
{
    corte_encontrado = true;
}
}

void CIgreedyB(vector<int> vetor_N1, vector<int> vetor_N2, vector<double> aux_denso, double &soma_demanda_cobertura, double B, vector<double> d, vector<int> &vetor_cobertura, bool &corte_encontrado, vector<double> x, double y_til,  vector<int> &best_CI, double &valor_best_CI)
{
vector<int> vetorDN1 = vetor_cobertura;
for (int i = 0; i < vetor_N1.size(); i++)
{
         vetorDN1.push_back(vetor_N1[i]);
}
vetor_cobertura.clear();
soma_demanda_cobertura = 0;
        int n = x.size();
        vector<double> aux2(n);
        for (int i = 0; i < n; i++)
        {
            aux2[i] = (y_til-x[i])/d[i];
        }
        funcao_ordena_crescente(vetorDN1, aux2);
        for (int i = 0; i < vetorDN1.size(); i++)
    {

    if (soma_demanda_cobertura <=  B  )
   {
       vetor_cobertura.push_back(vetorDN1[i]);
       soma_demanda_cobertura = soma_demanda_cobertura + d[vetorDN1[i]];
   }
   else
   {
       break;
   }
}   
  
if (soma_demanda_cobertura <= B)
{
funcao_ordena_decrescente(vetor_N2, d);
for (int i = 0; i < vetor_N2.size(); i++)
{
   if (soma_demanda_cobertura <=  B )
   {
       vetor_cobertura.push_back(vetor_N2[i]);

       soma_demanda_cobertura = soma_demanda_cobertura + d[vetor_N2[i]];
   }
   else
   {
       break;
   }
}
}
double soma_aux1 = 0;
for (int i = 0; i < vetor_cobertura.size(); i++)
{
     soma_aux1 = soma_aux1 + x[vetor_cobertura[i]];
}
double limitante_minimal = vetor_cobertura.size() - 1 ;
if (soma_aux1 - 0.000001 > limitante_minimal*y_til)
{
    corte_encontrado = true;
}

}

void AddOneCostumer(vector<int> vetor_N1, vector<int> vetor_N2, vector<double> aux_denso, double &soma_demanda_cobertura, double B, vector<double> d, vector<int> &vetor_cobertura, bool &corte_encontrado, vector<double> x ,double y_til,  vector<int> &best_CI, double &valor_best_CI  )
{

int melhor_cliente;
double x_melhor_cliente = 0;
int cliente_de_melhora;   //
if (y_til > 0.999999)
{
for (int i = 0; i < vetor_N1.size(); i++)
{
    if (  d[vetor_N1[i]] > (B - soma_demanda_cobertura) )  
    {
        cliente_de_melhora = vetor_N1[i];
        corte_encontrado = true;
        break;
    } 

}
if (corte_encontrado == true)
{
    vetor_cobertura.push_back(cliente_de_melhora);
    soma_demanda_cobertura = soma_demanda_cobertura + d[cliente_de_melhora]; 
}
}
else  // y_til < 1
{
double soma_lado_esquerdo = vetor_cobertura.size();
for (int i = 0; i < vetor_N1.size(); i++)
{
    if (  d[vetor_N1[i]] > (B - soma_demanda_cobertura)  && ( soma_lado_esquerdo + x[vetor_N1[i]] - 0.000001 > vetor_cobertura.size()*y_til ) )  
    {
        cliente_de_melhora = vetor_N1[i];
        corte_encontrado = true; //primeira melhora
        break;
     
    } 
}
if (corte_encontrado == true)
{
    vetor_cobertura.push_back(cliente_de_melhora);
    soma_demanda_cobertura = soma_demanda_cobertura + d[cliente_de_melhora]; 
}
}

} // fim void

void AddTwoCostumers(vector<int> vetor_N1, vector<int> vetor_N2, vector<double> aux_denso, double &soma_demanda_cobertura, double B, vector<double> d, vector<int> &vetor_cobertura, bool &corte_encontrado, vector<double> x, double y_til,  vector<int> &best_CI, double &valor_best_CI   )
{
int melhor_cliente1;
int melhor_cliente2;
double x_melhores_clientes = 0;
int cliente_de_melhora1;
int cliente_de_melhora2;
if (y_til > 0.999999)
{
for (int i1 = 0; i1 < vetor_N1.size() - 1; i1++)
{
    for (int i2 = i1 + 1; i2 < vetor_N1.size(); i2++)
    {
    if (  (d[vetor_N1[i1]] + d[vetor_N1[i2]] > (B - soma_demanda_cobertura)) && (   x[vetor_N1[i1]] + x[vetor_N1[i2]] - 1 > 0.000001    )   )
    {
        cliente_de_melhora1 = vetor_N1[i1];
        cliente_de_melhora2 = vetor_N1[i2];
        corte_encontrado = true;
        goto corte_encontrado_acrescenta2a; //first improvement
    }
    }
}

corte_encontrado_acrescenta2a:
if (corte_encontrado == true)
{
    vetor_cobertura.push_back(cliente_de_melhora1);
    vetor_cobertura.push_back(cliente_de_melhora2);
    soma_demanda_cobertura = soma_demanda_cobertura + d[cliente_de_melhora1] + d[cliente_de_melhora2]; 
}
}
else  //ytil < 1
{
double soma_lado_esquerdo = vetor_cobertura.size();
for (int i1 = 0; i1 < vetor_N1.size() - 1; i1++)
{
    for (int i2 = i1 + 1; i2 < vetor_N1.size(); i2++)
    {
    if (  (d[vetor_N1[i1]] + d[vetor_N1[i2]] > (B - soma_demanda_cobertura)) && (  soma_lado_esquerdo + x[vetor_N1[i1]] + x[vetor_N1[i2]] - 0.000001 > (vetor_cobertura.size() + 1)*y_til )   )
    {
        cliente_de_melhora1 = vetor_N1[i1];
        cliente_de_melhora2 = vetor_N1[i2];
        corte_encontrado = true;
        goto corte_encontrado_acrescenta2b;
    }
    }
}
corte_encontrado_acrescenta2b:
if (corte_encontrado == true)
{
    vetor_cobertura.push_back(cliente_de_melhora1);
    vetor_cobertura.push_back(cliente_de_melhora2);
    soma_demanda_cobertura = soma_demanda_cobertura + d[cliente_de_melhora1] + d[cliente_de_melhora2]; 
}

}
} //fim void

void AddThreeCostumers(vector<int> vetor_N1, vector<int> vetor_N2, vector<double> aux_denso, double &soma_demanda_cobertura, double B, vector<double> d, vector<int> &vetor_cobertura, bool &corte_encontrado, vector<double> x, double y_til,  vector<int> &best_CI, double &valor_best_CI   )
{
int melhor_cliente1;
int melhor_cliente2;
int melhor_cliente3;
double x_melhores_clientes = 0;
int cliente_de_melhora1;
int cliente_de_melhora2;
int cliente_de_melhora3;
if (y_til > 0.999999)
{
for (int i1 = 0; i1 < vetor_N1.size() - 1; i1++)
{
    for (int i2 = i1 + 1; i2 < vetor_N1.size(); i2++)
    {
    for (int i3 = i2 + 1; i3 < vetor_N1.size(); i3++)
    {
    if (  (d[vetor_N1[i1]] + d[vetor_N1[i2]] + d[vetor_N1[i3]] > (B - soma_demanda_cobertura)) && (   x[vetor_N1[i1]] + x[vetor_N1[i2]] + x[vetor_N1[i3]] - 2 > 0.000001    )   )
    {
        cliente_de_melhora1 = vetor_N1[i1];
        cliente_de_melhora2 = vetor_N1[i2];
        cliente_de_melhora3 = vetor_N1[i3];
        corte_encontrado = true;
        goto corte_encontrado_acrescenta3a; //first improvement
    }
    
    }
    }
}

corte_encontrado_acrescenta3a:
if (corte_encontrado == true)
{
    vetor_cobertura.push_back(cliente_de_melhora1);
    vetor_cobertura.push_back(cliente_de_melhora2);
    vetor_cobertura.push_back(cliente_de_melhora3);
    soma_demanda_cobertura = soma_demanda_cobertura + d[cliente_de_melhora1] + d[cliente_de_melhora2] + d[cliente_de_melhora3]; 
}
}
else  //ytil < 1
{
double soma_lado_esquerdo = vetor_cobertura.size();
for (int i1 = 0; i1 < vetor_N1.size() - 1; i1++)
{
    for (int i2 = i1 + 1; i2 < vetor_N1.size(); i2++)
    {
    for (int i3 = i2 + 1; i3 < vetor_N1.size(); i3++)
    {
    if (  (d[vetor_N1[i1]] + d[vetor_N1[i2]] + d[vetor_N1[i3]]  > (B - soma_demanda_cobertura)) && (  soma_lado_esquerdo + x[vetor_N1[i1]] + x[vetor_N1[i2]] + x[vetor_N1[i3]] - 0.000001 > (vetor_cobertura.size() + 2)*y_til )   )
    {
        cliente_de_melhora1 = vetor_N1[i1];
        cliente_de_melhora2 = vetor_N1[i2];
        cliente_de_melhora3 = vetor_N1[i3];
        corte_encontrado = true;
        goto corte_encontrado_acrescenta3b;
    }
    }
    }
}
corte_encontrado_acrescenta3b:

if (corte_encontrado == true)
{
    vetor_cobertura.push_back(cliente_de_melhora1);
    vetor_cobertura.push_back(cliente_de_melhora2);
    vetor_cobertura.push_back(cliente_de_melhora3);
    soma_demanda_cobertura = soma_demanda_cobertura + d[cliente_de_melhora1] + d[cliente_de_melhora2] + d[cliente_de_melhora3];
}
}
}

void FCColumnGeneration (vector<double> Pesos , vector<int> Itens, double B, vector<double> &coef, vector<int> &clientes_corte, double &limitante, vector<double> &x, bool &corte_violado, bool &impossivel_cortar, int &tentativas_opcao0, int &total_tentativas, int bigM, int max_iter_fc)
{
///Reduz o problema da mochila; 
vector<int> J1;
vector <int> J0;
vector<int> Jf;
vector<double> xf;
vector<double> Pesos_correspondentes;
bool duais_fracionarios = false;
for (int i = 0; i < Itens.size(); i++)
{
    if (x[Itens[i]] > 0.99999)
    {
        J1.push_back(Itens[i]);
         B = B - Pesos[Itens[i]] ;
    }
  
}
for (int i = 0; i < Itens.size(); i++)
{
    if (x[Itens[i]] > 0.99999) continue;

     
            if ( (x[Itens[i]] < 0.00001) )
            {
            J0.push_back(Itens[i]);
             }
             else
                {
                   Jf.push_back(Itens[i]);
                  Pesos_correspondentes.push_back(Pesos[Itens[i]]);
                  xf.push_back(x[Itens[i]]);
                }
    }    
vector<int> itens_mochila(Jf.size());
    for (int i = 0; i < Jf.size(); i++)
    {
       itens_mochila[i] = i;
    }

IloEnv env;
IloModel mod(env);
IloCplex cplex(mod);
IloObjective Objetivo = IloAdd(mod, IloMinimize(env));
IloNumArray lado_direito(env,Jf.size());

for (int i = 0; i < Jf.size(); i++)
{
lado_direito[i] = x[Jf[i]];     
}
IloRangeArray Restricoes = IloAdd(mod, IloRangeArray(env, lado_direito, IloInfinity));
IloNumVarArray Variaveis(env);
IloNumColumnArray Colunas(env);
for (IloInt j = 0; j < Restricoes.getSize(); j++)
{
IloNumColumn col = Objetivo(1);
for (IloInt i = 0; i < Restricoes.getSize(); ++i)
{
    if(i == j)
    {
     col += Restricoes[i](1);
    }
}
IloNumVar var(col, 0,  IloInfinity);
Variaveis.add(var);
Colunas.add(col);

}
    cplex.setWarning(env.getNullStream()); // Eliminar warnings
    cplex.setOut(env.getNullStream()); //Eliminar os logs do solver
    cplex.setParam(IloCplex::Threads,1);
    cplex.setParam(IloCplex::Param::Parallel,1);
cplex.solve();
vector<double> duais(Restricoes.getSize());
vector<double> vetor_duais_fracionarios;
for ( int i = 0; i < Restricoes.getSize(); i++)
{
     duais[i] = cplex.getDual(Restricoes[i]) ;
     if ( abs(round(duais[i]) - duais[i]) >= 0.000001) 
     {
     vetor_duais_fracionarios.push_back(duais[i]);
     }
     else 
     {
     duais[i] = round(duais[i]);      
     } 
       
}
if (vetor_duais_fracionarios.size() > 0.5) 
{
duais_fracionarios = true;
}
vector<int> itens_otimo;  // itens da solucao otima do problema da mochila indica as posicoes não nulas da coluna a ser adicionada no problema
double otimo_mochila;
double otimo_mochila_antes = 0; 
vector<double> duais_antes = duais;   //guarda o duais
vector<double> duais_copia = duais;
int conta_iter = 0;
while (  (corte_violado == false)  )
{
conta_iter = conta_iter + 1;
if ( (conta_iter > max_iter_fc) && (cplex.getObjValue() > 1.1 ))
{
    impossivel_cortar = true;
    return;
}
if ( cplex.getObjValue() <= 1 )
{
    impossivel_cortar = true;
    return;
}
else //o corte é violado
{
bool nao_alterou_corte = false; 
duais_antes = duais;
duais_copia = duais; 
limitante = 1;  
    if (duais_fracionarios == true) //duais sao fracionarios entao ha duas opcoes de arredondamento apos multiplicacao por bigM
    {
    total_tentativas = total_tentativas + 1;    
    bool opcao0 = true; // kaparis2010
    bool opcao1 = false; // gadegaard2018
    bool opcao2 = false; // yang2012
    bool opcao3 = false; // problema inteiro boccia
    bool opcao4 = false; // nao fazer arredondamento. 
    nao_alterou_corte = false;
    bool opcao0_violada = false;
    bool opcao1_violada = false;
    bool opcao2_violada = false;
    bool opcao3_violada = false;
    bool opcao4_violada = false;
    bool opcao0_valida = false;
    bool opcao1_valida = false;
    bool opcao2_valida = false;
    bool opcao3_valida = false;
    bool opcao4_valida = false;
    if (opcao0 == true)
    {
    opcao0_violada = true;
    double menor_valor_fracionario = 100000000;
    for (int i = 0; i < vetor_duais_fracionarios.size(); i++)
    {
    if (vetor_duais_fracionarios[i] < menor_valor_fracionario)
    {
            menor_valor_fracionario = vetor_duais_fracionarios[i];
    }
    }
    for (int i = 0; i < duais_copia.size(); i++)
    {
    if ( abs(round(duais_copia[i]/menor_valor_fracionario) - duais_copia[i]/menor_valor_fracionario ) < 0.000001) 
    {
        duais_copia[i] = round(duais_copia[i]/menor_valor_fracionario);
    }
    else 
    {
        opcao0_violada = false;
        break;    
    }
    }
if (opcao0_violada == true)
{
    if (abs(round((double) 1/menor_valor_fracionario) - 1/menor_valor_fracionario  ) >= 0.000001    )
    {
        opcao0_violada = false;
    }
    else
    {
        limitante = round((double) 1/menor_valor_fracionario); 
    }
    }
    } // fim opcao 0 
    if (opcao0_violada == true )  // verifica se o corte é válido 
    {
       tentativas_opcao0 = tentativas_opcao0 + 1;
       duais = duais_copia;
    vector<int> itens_otimo_copia = itens_otimo;    
    itens_otimo.clear();
    vector<int> indices_considerados; // itens mochila com duais (coeficientes nao nulos)
    vector<double> Pesos_correspondentes2; // pesos correspondentes dos itens da mochila com coeficientes nao nulos
    vector<double> duais2; // duais correspondentes dos itens da mochila com coeficientes nao nulos
    vector<int> itens_mochila2;
    vector<int> indices_otimo;
    double soma_indices_considerados = 0;
    double soma_coeficientes_fo = 0;
    for (int i = 0; i < duais.size(); i++)
    {
        soma_coeficientes_fo = soma_coeficientes_fo + duais[i];
    }
    

     for (int i = 0; i < itens_mochila.size(); i++)
    {
        if (duais[i] > 0 ) 
        {
            indices_considerados.push_back(i);
            Pesos_correspondentes2.push_back(Pesos_correspondentes[i]);
            duais2.push_back(duais[i]);
            soma_indices_considerados = soma_indices_considerados + Pesos_correspondentes[i]; 
        }
    }
    for (int i = 0; i < indices_considerados.size(); i++)
    {
        itens_mochila2.push_back(i);
    }
                        ///////
        if ( (soma_indices_considerados <= B) /*&& (4 > 5)*/)
        {
            otimo_mochila = soma_coeficientes_fo;
            for (int i = 0; i < indices_considerados.size(); i++)
            {
                itens_otimo.push_back(itens_mochila[indices_considerados[i]]);
            }
        }
        else
        {
         SolveKnapsack( B, duais2, itens_mochila2 , Pesos_correspondentes2, otimo_mochila, indices_otimo) ;
        itens_otimo.clear();
        for (int i = 0; i < indices_otimo.size(); i++)
            {
                itens_otimo.push_back(indices_considerados[indices_otimo[i]]);
            }
        }
    if (otimo_mochila <= limitante)
    {
      corte_violado = true;
      clientes_corte = Jf;
      coef = duais;
      limitante = otimo_mochila;
    }
    else //adiciona coluna se alterou o corte
    {
            if (itens_otimo == itens_otimo_copia)
        {
            nao_alterou_corte = true;
        }
    if (nao_alterou_corte == false) 
    { 
    IloNumColumn col = Objetivo(1);
    for (int i = 0; i < itens_otimo.size(); i++)
        {
            col += Restricoes[itens_otimo[i]](1);
        }
    IloNumVar var(col, 0,  IloInfinity);
    Variaveis.add(var);
    Colunas.add(col);
        cplex.setWarning(env.getNullStream()); // Eliminar warnings
    cplex.setOut(env.getNullStream()); //Eliminar os logs do solver
    cplex.setParam(IloCplex::Threads,1);
    cplex.setParam(IloCplex::Param::Parallel,1);
    cplex.solve();
     duais_fracionarios = false;
    vetor_duais_fracionarios.clear();
    for ( int i = 0; i < Restricoes.getSize(); i++)
    {
     duais[i] = cplex.getDual(Restricoes[i]) ;
     if ( abs(round(duais[i]) - duais[i]) >= 0.000001) 
     {
     vetor_duais_fracionarios.push_back(duais[i]);
     }
     else 
     {
     duais[i] = round(duais[i]);      
     }   
    }
    if (vetor_duais_fracionarios.size() > 0.5) 
    {
    duais_fracionarios = true;
    }
    }
    else //nao alterou o corte 
    {
          opcao1 = true;
    }
    }
    }
    else
    {
        opcao1 = true;
    }
    if (opcao1 == true)
    {
    duais = duais_antes;
    duais_copia = duais_antes;
    opcao1_violada = true;
    nao_alterou_corte = false;    
    for (int i = 0; i < Jf.size(); i++)
    {
        duais_copia[i] = round(bigM*duais_copia[i]);
    }
    limitante = bigM; 
    double aux = 0;
    for (int i = 0; i < Jf.size(); i++)
    {
    aux = aux + duais_copia[i]*xf[i];
    }
    if (aux > limitante)
    {
     opcao1_violada = true;
    }
    else
    {
        opcao1_violada = false;
    }
    }
    if ((opcao1 == true) && (opcao1_violada == true)  )
    {
    duais = duais_copia;
    vector<int> itens_otimo_copia = itens_otimo;    
    itens_otimo.clear();
    vector<int> indices_considerados; // itens mochila com duais (coeficientes nao nulos)
    vector<double> Pesos_correspondentes2; // pesos correspondentes dos itens da mochila com coeficientes nao nulos
    vector<double> duais2; // duais correspondentes dos itens da mochila com coeficientes nao nulos
    vector<int> itens_mochila2;
    vector<int> indices_otimo;
    double soma_indices_considerados = 0;
    double soma_coeficientes_fo = 0;
    for (int i = 0; i < duais.size(); i++)
    {
        soma_coeficientes_fo = soma_coeficientes_fo + duais[i];
    }
     for (int i = 0; i < itens_mochila.size(); i++)
    {
        if (duais[i] > 0 ) 
        {
            indices_considerados.push_back(i);
            Pesos_correspondentes2.push_back(Pesos_correspondentes[i]);
            duais2.push_back(duais[i]);
            soma_indices_considerados = soma_indices_considerados + Pesos_correspondentes[i]; 
        }
    }
    for (int i = 0; i < indices_considerados.size(); i++)
    {
        itens_mochila2.push_back(i);
    }
        if (soma_indices_considerados <= B)
        {
            otimo_mochila = soma_coeficientes_fo;
            for (int i = 0; i < indices_considerados.size(); i++)
            {
                itens_otimo.push_back(itens_mochila[indices_considerados[i]]);
            }
        }
        else
        {
         SolveKnapsack( B, duais2, itens_mochila2 , Pesos_correspondentes2, otimo_mochila, indices_otimo) ;
        itens_otimo.clear();
        for (int i = 0; i < indices_otimo.size(); i++)
            {
                itens_otimo.push_back(indices_considerados[indices_otimo[i]]);
            }
        }
                                    if (itens_otimo == itens_otimo_copia)
        {
            nao_alterou_corte = true;
        }
    if (otimo_mochila <= limitante)
    {
      corte_violado = true;
      clientes_corte = Jf;
      coef = duais;
        limitante = otimo_mochila;
    }
    else //adiciona coluna só adiciona coluna se tiver alterado o corte 
    {
        if ( (otimo_mochila - limitante) < 0.00001*limitante )
        {
         opcao2 =true;
        }   
                            if ( (itens_otimo == itens_otimo_copia) && (corte_violado == false))
        {
        nao_alterou_corte = true;
        }
    if ( (nao_alterou_corte == false) && (opcao2 == false))
    {
    IloNumColumn col = Objetivo(1);
    for (int i = 0; i < itens_otimo.size(); i++)
        {
            col += Restricoes[itens_otimo[i]](1);
        }
        
    IloNumVar var(col, 0,  IloInfinity);
    Variaveis.add(var);
    Colunas.add(col);
        cplex.setWarning(env.getNullStream()); // Eliminar warnings
    cplex.setOut(env.getNullStream()); //Eliminar os logs do solver
    cplex.setParam(IloCplex::Threads,1);
    cplex.setParam(IloCplex::Param::Parallel,1);
    cplex.solve();
     duais_fracionarios = false;
    vetor_duais_fracionarios.clear();
    for ( int i = 0; i < Restricoes.getSize(); i++)
    {
     duais[i] = cplex.getDual(Restricoes[i]) ;
     if ( abs(round(duais[i]) - duais[i]) >= 0.000001) 
     {
     vetor_duais_fracionarios.push_back(duais[i]);
     }
     else 
     {
     duais[i] = round(duais[i]);      
     }   
    }
    if (vetor_duais_fracionarios.size() > 0.5) 
    {
    duais_fracionarios = true;
    }  
    }
    else  // nao alterou o corte
    {
       opcao2 = true;
    }

    }
    } // fim opcao1 deu certo 
if ( (opcao1 == true) && (opcao1_violada == false) )
{
    opcao2 = true;
}
    if (opcao2 == true)
    {
        opcao2_violada = true;
        duais = duais_antes;
        duais_copia = duais_antes;
    nao_alterou_corte = false;    
            for (int i = 0; i < Jf.size(); i++)
    {
        duais_copia[i] = duais_copia[i]*bigM;
        if (abs(duais_copia[i] - round(duais_copia[i])) < 0.000001) duais_copia[i] = round(duais_copia[i]);
        duais_copia[i] = floor(duais_copia[i]);
    }
    limitante = bigM; 
    double aux = 0;
    for (int i = 0; i < Jf.size(); i++)
    {
    aux = aux + duais_copia[i]*xf[i];
    }
    if (aux > limitante)
    {
        opcao2_violada = true;
    }
    else
    {
        opcao2_violada = false;
        impossivel_cortar = true;
        return;
    }
    }
    if  ( (opcao2 == true) && (opcao2_violada == true) )
    {
        duais = duais_copia;
    vector<int> itens_otimo_copia = itens_otimo;    
    itens_otimo.clear();
    vector<int> indices_considerados; // itens mochila com duais (coeficientes nao nulos)
    vector<double> Pesos_correspondentes2; // pesos correspondentes dos itens da mochila com coeficientes nao nulos
    vector<double> duais2; // duais correspondentes dos itens da mochila com coeficientes nao nulos
    vector<int> itens_mochila2;
    vector<int> indices_otimo;
    double soma_indices_considerados = 0;
    double soma_coeficientes_fo = 0;
    for (int i = 0; i < duais.size(); i++)
    {
        soma_coeficientes_fo = soma_coeficientes_fo + duais[i];
    }
     for (int i = 0; i < itens_mochila.size(); i++)
    {
        if (duais[i] > 0 ) 
        {
            indices_considerados.push_back(i);
            Pesos_correspondentes2.push_back(Pesos_correspondentes[i]);
            duais2.push_back(duais[i]);
            soma_indices_considerados = soma_indices_considerados + Pesos_correspondentes[i]; 
        }
    }
    for (int i = 0; i < indices_considerados.size(); i++)
    {
        itens_mochila2.push_back(i);
    }
        if (soma_indices_considerados <= B)
        {
            otimo_mochila = soma_coeficientes_fo;
            for (int i = 0; i < indices_considerados.size(); i++)
            {
                itens_otimo.push_back(itens_mochila[indices_considerados[i]]);
            }
                    if (itens_otimo == itens_otimo_copia)
        {
            nao_alterou_corte = true;
            impossivel_cortar = true;
            return;
        }
        }
        else
        {
         SolveKnapsack( B, duais2, itens_mochila2 , Pesos_correspondentes2, otimo_mochila, indices_otimo) ;
        itens_otimo.clear();
        for (int i = 0; i < indices_otimo.size(); i++)
            {
                itens_otimo.push_back(indices_considerados[indices_otimo[i]]);
            }
                    if (itens_otimo == itens_otimo_copia)
        {
            nao_alterou_corte = true;
            impossivel_cortar = true;
            return;
        }
        }
    if (otimo_mochila <= limitante)
    {
      corte_violado = true;
      clientes_corte = Jf;
      coef = duais;
        limitante = otimo_mochila;
    }
    else //adiciona coluna 
    {
    IloNumColumn col = Objetivo(1);
    for (int i = 0; i < itens_otimo.size(); i++)
        {
            col += Restricoes[itens_otimo[i]](1);
        }
    IloNumVar var(col, 0,  IloInfinity);
    Variaveis.add(var);
    Colunas.add(col);
        cplex.setWarning(env.getNullStream()); // Eliminar warnings
    cplex.setOut(env.getNullStream()); //Eliminar os logs do solver
    cplex.setParam(IloCplex::Threads,1);
    cplex.setParam(IloCplex::Param::Parallel,1);
    cplex.solve();
     duais_fracionarios = false;
    vetor_duais_fracionarios.clear();
    for ( int i = 0; i < Restricoes.getSize(); i++)
    {
     duais[i] = cplex.getDual(Restricoes[i]) ;
     if ( abs(round(duais[i]) - duais[i]) >= 0.000001) 
     {
     vetor_duais_fracionarios.push_back(duais[i]);
     }
     else 
     {
     duais[i] = round(duais[i]);      
     }   
    }
    if (vetor_duais_fracionarios.size() > 0.5) 
    {
    duais_fracionarios = true;
    }  
    }
    } // fim opcao 2 deu certo
    }
   else  //duais sao inteiros entao nenhum arredondamento é necessário. 
    {
                            itens_otimo.clear();
                            vector<int> itens_otimo_copia = itens_otimo;
    vector<int> indices_considerados; // itens mochila com duais (coeficientes nao nulos)
    vector<double> Pesos_correspondentes2; // pesos correspondentes dos itens da mochila com coeficientes nao nulos
    vector<double> duais2; // duais correspondentes dos itens da mochila com coeficientes nao nulos
    vector<int> itens_mochila2;
    vector<int> indices_otimo;
    double soma_indices_considerados = 0;
    double soma_coeficientes_fo = 0;
    for (int i = 0; i < duais.size(); i++)
    {
        soma_coeficientes_fo = soma_coeficientes_fo + duais[i];
    }
     for (int i = 0; i < itens_mochila.size(); i++)
    {
        if (duais[i] > 0 ) 
        {
            indices_considerados.push_back(i);
            Pesos_correspondentes2.push_back(Pesos_correspondentes[i]);
            duais2.push_back(duais[i]);
            soma_indices_considerados = soma_indices_considerados + Pesos_correspondentes[i]; 
        }
    }
    for (int i = 0; i < indices_considerados.size(); i++)
    {
        itens_mochila2.push_back(i);
    }
        if (soma_indices_considerados <= B)
        {
            otimo_mochila = soma_coeficientes_fo;
            for (int i = 0; i < indices_considerados.size(); i++)
            {
                itens_otimo.push_back(itens_mochila[indices_considerados[i]]);
            }
              if (itens_otimo == itens_otimo_copia)
        {
            nao_alterou_corte = true;
        }
        }
        else
        {
         SolveKnapsack( B, duais2, itens_mochila2 , Pesos_correspondentes2, otimo_mochila, indices_otimo) ;
        itens_otimo.clear();
        for (int i = 0; i < indices_otimo.size(); i++)
            {
                itens_otimo.push_back(indices_considerados[indices_otimo[i]]);
            }
              if (itens_otimo == itens_otimo_copia)
        {
            nao_alterou_corte = true;
        }
        }
       if (nao_alterou_corte == false)
       { 
    if (otimo_mochila <= 1)
    {
      corte_violado = true;
      clientes_corte = Jf;
      coef = duais;
      limitante = 1;
    }
    else //adiciona coluna 
    {
    IloNumColumn col = Objetivo(1);
    for (int i = 0; i < itens_otimo.size(); i++)
        {
            col += Restricoes[itens_otimo[i]](1);
        }
    IloNumVar var(col, 0,  IloInfinity);
    Variaveis.add(var);
    Colunas.add(col);
        cplex.setWarning(env.getNullStream()); // Eliminar warnings
    cplex.setOut(env.getNullStream()); //Eliminar os logs do solver
    cplex.setParam(IloCplex::Threads,1);
    cplex.setParam(IloCplex::Param::Parallel,1);
    cplex.solve();
       duais_fracionarios = false;
    vetor_duais_fracionarios.clear();
    for ( int i = 0; i < Restricoes.getSize(); i++)
    {
     duais[i] = cplex.getDual(Restricoes[i]) ;
     if ( abs(round(duais[i]) - duais[i]) >= 0.000001) 
     {
     vetor_duais_fracionarios.push_back(duais[i]);
     }
     else 
     {
     duais[i] = round(duais[i]);      
     }   
    }
    if (vetor_duais_fracionarios.size() > 0.5) 
    {
    duais_fracionarios = true;
    }
    }
    }
    } //fim duais nao sao fracionarios 
}
} // fim loop while fenchel cut 
if ( corte_violado == true)   
{
if (duais_fracionarios == false) //lifting é feito com limitante 1
{
if (  (J1.size() > 0)   ) 
{
for (int i = 0; i < J1.size(); i++)
{
vector<int> itens_mochila;
vector<double> pesos_mochila;
vector<double> valores_mochila;
for (int i2 = 0; i2 < clientes_corte.size(); i2++)
{
    itens_mochila.push_back(i2);
    pesos_mochila.push_back(Pesos[clientes_corte[i2]]);
    //valores_mochila.push_back(coef[i]);
}
                for (int i2 = 0; i2 < coef.size(); i2++)
        {
                    valores_mochila.push_back(coef[i2]);
        }
B  = B +  Pesos[J1[i]];
vector<int> itens_otimo;
double otimo_mochila = 0;
SolveKnapsack( B, valores_mochila, itens_mochila , pesos_mochila, otimo_mochila, itens_otimo) ;
clientes_corte.push_back(J1[i]);
double coef_aux = otimo_mochila - limitante;
coef.push_back(coef_aux); 
limitante = limitante + coef_aux;
}
}
if ( (J0.size() > 0)   )  
{
bool opcao1 = false; //faz up lifting utilizando a estrategia ECI 
if (opcao1 == true)
{
int maior_peso = 0;
int indice_cliente_de_maior_peso; 
for (int i = 0; i < clientes_corte.size(); i++)
{
    if (Pesos[clientes_corte[i]] > maior_peso )
    {
        maior_peso = Pesos[clientes_corte[i]]; 
        indice_cliente_de_maior_peso = i;
    }
}
for (int i = 0; i < J0.size(); i++)
{
    if (Pesos[J0[i]]>= maior_peso)
    {
        clientes_corte.push_back(J0[i]);
        coef.push_back(coef[indice_cliente_de_maior_peso]);
    }
}
}
else
{
for (int i = 0; i < J0.size(); i++)
{
vector<int> itens_mochila;
vector<double> pesos_mochila;
vector<double> valores_mochila;
for (int i2 = 0; i2 < clientes_corte.size(); i2++)
{
    itens_mochila.push_back(i2);
    pesos_mochila.push_back(Pesos[clientes_corte[i2]]);
    valores_mochila.push_back(coef[i2]);
}
B = B - Pesos[J0[i]]; // qual é o valor máximo considerando as variáveis do corte caso a variável presente entre no corte com valor 1.  
if(B < 0) continue;
vector<int> itens_otimo;
double otimo_mochila = 0;
SolveKnapsack( B, valores_mochila, itens_mochila , pesos_mochila, otimo_mochila, itens_otimo) ;
int coef_aux = round(limitante - otimo_mochila);
if (coef_aux != 0)
{
coef.push_back(coef_aux); 
clientes_corte.push_back(J0[i]);
}
B = B +  Pesos[J0[i]];   //volta a capacidade que estava antes
}
}
} // fim j0.size > 0
}
else  //lifting é feito com limitante bigM
{
if (  (J1.size() > 0)   ) 
{
for (int i = 0; i < J1.size(); i++)
{
vector<int> itens_mochila;
vector<double> pesos_mochila;
vector<double> valores_mochila;
for (int i2 = 0; i2 < clientes_corte.size(); i2++)
{
    itens_mochila.push_back(i2);
    pesos_mochila.push_back(Pesos[clientes_corte[i2]]);
    //valores_mochila.push_back(coef[i]);
}
                for (int i2 = 0; i2 < coef.size(); i2++)
        {
                    valores_mochila.push_back(coef[i2]);
        }
B  = B +  Pesos[J1[i]];
vector<int> itens_otimo;
double otimo_mochila = 0;
SolveKnapsack( B, valores_mochila, itens_mochila , pesos_mochila, otimo_mochila, itens_otimo) ;
clientes_corte.push_back(J1[i]);
double coef_aux = otimo_mochila - limitante;
coef.push_back(coef_aux); 
limitante = limitante + coef_aux;
}
}
if ( (J0.size() > 0)   )  
{
bool opcao1 = false;   // utiliza estrategia ECI 
if (opcao1 == true)
{
int maior_peso = 0;
int indice_cliente_de_maior_peso; 

for (int i = 0; i < clientes_corte.size(); i++)
{
    if (Pesos[clientes_corte[i]] > maior_peso )
    {
        maior_peso = Pesos[clientes_corte[i]]; 
        indice_cliente_de_maior_peso = i;
    }
}
for (int i = 0; i < J0.size(); i++)
{
    if (Pesos[J0[i]]>= maior_peso)
    {
        clientes_corte.push_back(J0[i]);
        coef.push_back(coef[indice_cliente_de_maior_peso]);
    }
}
}
else
{
for (int i = 0; i < J0.size(); i++)
{
vector<int> itens_mochila;
vector<double> pesos_mochila;
vector<double> valores_mochila;
for (int i2 = 0; i2 < clientes_corte.size(); i2++)
{
    itens_mochila.push_back(i2);
    pesos_mochila.push_back(Pesos[clientes_corte[i2]]);
    valores_mochila.push_back(coef[i2]);
}
B = B - Pesos[J0[i]]; // qual é o valor máximo considerando as variáveis do corte caso a variável presente entre no corte com valor 1.  
if(B < 0) continue;
vector<int> itens_otimo;
double otimo_mochila = 0;
SolveKnapsack( B, valores_mochila, itens_mochila , pesos_mochila, otimo_mochila, itens_otimo) ;
int coef_aux = round(limitante - otimo_mochila);
if (coef_aux != 0)
{
coef.push_back(coef_aux); 
clientes_corte.push_back(J0[i]);
}
B = B +  Pesos[J0[i]];   //volta a capacidade que estava antes
}
} // fim else
}  //fim J0.size() > 0
} // fim void aplica fenchel cut novo
}
}
void CIexact(double B, vector<int> vetor_N1, vector<int> vetor_N2, vector<int> vetor_D , double y_til, vector<double> x, vector<double> d, vector<int> &vetor_cobertura_auxiliar, double &soma_demandas_cobertura_auxiliar, bool &corte_encontrado, double soma_demandas_D, bool &cobertura_encontrada,  vector<int> &best_CI, double &valor_best_CI)
{
if (y_til > 0.999999)
{
int k_best = -1;
int inf = pow(10,10);
double g_best = inf;
    list<int> itens_adicioandos;
B = B - soma_demandas_D;
int n = vetor_N1.size();
vector<vector<double>> f(n+1, vector<double>(B+1));  
vector<double> g(n+1);
for (int k = 0; k < n+1; k++)
{
    g[k] = inf;
}
for (int r = 0; r < B+1 ; r++)
{
    f[0][r] = inf;
}
f[0][0] = 0;
for (int k = 1; k < n+1 ; k++)
{
    for (int r = 0; r < B+1 ; r++)
    {
        f[k][r] = inf;
    }
    
}

for (int k = 1; k < n+1; k++)
{
    for (int r = 0; r < B+1; r++)
    {
        if ( f[k-1][r] < f[k][r])
        {
            f[k][r] = f[k-1][r];
        }
    }
    for (int r = 0; r < B - d[vetor_N1[k-1]] + 1; r++)
    {
        if ( f[k-1][r] + (y_til - x[vetor_N1[k-1]] ) < f[k][r + d[vetor_N1[k-1] ]])
        {
            f[k][r + d[vetor_N1[k-1]]] = f[k-1][r] + (y_til - x[vetor_N1[k-1]] );
        }
    }
        for (int r = B - d[vetor_N1[k-1]] + 1 ; r < B + 1 ; r++)
    {
        if ( f[k-1][r] + (y_til -x[vetor_N1[k-1]] ) < g[k] )
        {
            g[k] = f[k-1][r] + (y_til - x[vetor_N1[k-1]] );
            if (g[k] < g_best)
            {
                k_best = k;
                g_best = g[k];
            }   
        }
    }
    if (g[k] < y_til - 0.000001)
    {
        itens_adicioandos.push_back(vetor_N1[k-1]);
        soma_demandas_cobertura_auxiliar = soma_demandas_cobertura_auxiliar + d[vetor_N1[k - 1]];
        double demanda_minima = B + 1 - d[vetor_N1[k-1]];
        double soma_x = g[k] - (1 -x[vetor_N1[k-1]] ); 
       if (  (soma_x  < 0.000001) && (demanda_minima < 0) ) 
        {
                       itens_adicioandos.reverse();
        for (std::list<int>::iterator k1 = itens_adicioandos.begin(); k1 != itens_adicioandos.end(); k1++ ){
        vetor_cobertura_auxiliar.push_back(*k1);
        }
            corte_encontrado = true;
            cobertura_encontrada = true;
            return;
        }
        int x_pivo = k-1;
        int y_pivo;
        for (int r = demanda_minima; r < B+1; r++)
        {
           if (  (f[x_pivo][r] - soma_x) < 0.000001  ) 
           {
               y_pivo = r;
               break;
           } 
        }
        while (x_pivo >=1 )
        {
         if (x_pivo > 0)
         {
         if (  (f[x_pivo - 1][y_pivo] == f[x_pivo][y_pivo])  )
         {
            x_pivo = x_pivo - 1;
         }
         else
         {
             itens_adicioandos.push_back(vetor_N1[x_pivo-1]);
              soma_demandas_cobertura_auxiliar = soma_demandas_cobertura_auxiliar + d[vetor_N1[x_pivo - 1]];
              demanda_minima = demanda_minima - d[vetor_N1[x_pivo-1]];
              soma_x = soma_x - ( y_til - x[vetor_N1[x_pivo-1]]  ) ;
              if (  (soma_x  < 0.000001) && (demanda_minima < 0) ) 
            {
          itens_adicioandos.reverse();
        for (std::list<int>::iterator k1 = itens_adicioandos.begin(); k1 != itens_adicioandos.end(); k1++ ){
        vetor_cobertura_auxiliar.push_back(*k1);
        }
            corte_encontrado = true;
            cobertura_encontrada = true;
            return;
            }
              x_pivo = x_pivo - 1;
                for (int r = demanda_minima; r < B+1; r++)
                {
                if (abs (f[x_pivo][r] - soma_x) < 0.000001 ) 
                {
                    y_pivo = r;
                break;
                } 
                }
         }
         }
        }// fim while 
        itens_adicioandos.reverse();
        for (std::list<int>::iterator k = itens_adicioandos.begin(); k != itens_adicioandos.end(); k++ ){
        vetor_cobertura_auxiliar.push_back(*k);
        }
        corte_encontrado = true;
        cobertura_encontrada = true;
        return;   
    } // termina if o corte é encontrado
}   // fim for k = 1 ...
if (g_best < inf)
{
    cobertura_encontrada = true;
}
if ( (corte_encontrado == false) && (cobertura_encontrada == true)   )
{
itens_adicioandos.push_back(vetor_N1[k_best-1]);
        soma_demandas_cobertura_auxiliar = soma_demandas_cobertura_auxiliar + d[vetor_N1[k_best - 1]];
        double demanda_minima = B + 1 - d[vetor_N1[k_best-1]];
        double soma_x = g_best - (1 -x[vetor_N1[k_best-1]] ); 
       if (  (soma_x  < 0.000001) && (demanda_minima < 0) ) 
        {
                       itens_adicioandos.reverse();
        for (std::list<int>::iterator k1 = itens_adicioandos.begin(); k1 != itens_adicioandos.end(); k1++ ){
        vetor_cobertura_auxiliar.push_back(*k1);
        }
            return;
        }
        int x_pivo = k_best-1;
        int y_pivo;
        for (int r = demanda_minima; r < B+1; r++)
        {
           if (  (f[x_pivo][r] - soma_x) < 0.000001  ) 
           {
               y_pivo = r;
               break;
           } 
        }
        while (x_pivo >=1 )
        {
         if (x_pivo > 0)
         {
         if (  (f[x_pivo - 1][y_pivo] == f[x_pivo][y_pivo])  )
         {
            x_pivo = x_pivo - 1;
         }
         else
         {
             itens_adicioandos.push_back(vetor_N1[x_pivo-1]);
              soma_demandas_cobertura_auxiliar = soma_demandas_cobertura_auxiliar + d[vetor_N1[x_pivo - 1]];
              demanda_minima = demanda_minima - d[vetor_N1[x_pivo-1]];
              soma_x = soma_x - ( y_til - x[vetor_N1[x_pivo-1]]  ) ;
              if (  (soma_x  < 0.000001) && (demanda_minima < 0) ) 
            {
                           itens_adicioandos.reverse();
        for (std::list<int>::iterator k1 = itens_adicioandos.begin(); k1 != itens_adicioandos.end(); k1++ ){
        vetor_cobertura_auxiliar.push_back(*k1);
        }
            return;
            }
              x_pivo = x_pivo - 1;
                for (int r = demanda_minima; r < B+1; r++)
                {
                if (abs (f[x_pivo][r] - soma_x) < 0.000001 ) 
                {
                    y_pivo = r;
                break;
                } 
                }
         }
         }
        }// fim while 
        itens_adicioandos.reverse();
        for (std::list<int>::iterator k = itens_adicioandos.begin(); k != itens_adicioandos.end(); k++ ){
        vetor_cobertura_auxiliar.push_back(*k);
        }
        return;
} /// fim corte encontrado é false
}  /// fim y_til = 1
else
{
int inf = pow(10,10);    
int k_best = -1;
double g_best = inf;    
list<int> itens_adicioandos;
vetor_cobertura_auxiliar.clear();
soma_demandas_cobertura_auxiliar = 0;
funcao_ordena_decrescente(vetor_N2,d);
funcao_ordena_decrescente2(vetor_N1, x, d);
vector<int> itens = vetor_D;
for (int i = 0; i < vetor_N1.size(); i++)
{
    itens.push_back(vetor_N1[i]);
}
for (int i = 0; i < vetor_N2.size(); i++)
{
    if (d[vetor_N2[i]] <= B) itens.push_back(vetor_N2[i]);
}
int n = itens.size();
vector<vector<double>> f(n+1, vector<double>(B+1));  
vector<double> g(n+1);
for (int k = 0; k < n+1; k++)
{
    g[k] = inf;
}
for (int r = 0; r < B+1 ; r++)
{
    f[0][r] = inf;
}
f[0][0] = 0;
for (int k = 1; k < n+1 ; k++)
{
    for (int r = 0; r < B+1 ; r++)
    {
        f[k][r] = inf;
    }
}
for (int k = 1; k < n+1; k++)
{
    for (int r = 0; r < B+1; r++)
    {
        if ( f[k-1][r] < f[k][r])
        {
            f[k][r] = f[k-1][r];
        }
    }
    for (int r = 0; r < B - d[itens[k-1]] + 1; r++)
    {
        if ( f[k-1][r] + (y_til -x[itens[k-1]] ) < f[k][r + d[itens[k-1]]])
        {
            f[k][r + d[itens[k-1]]] = f[k-1][r] + (y_til -x[itens[k-1] ]);
        }
    }
        for (int r = B - d[itens[k-1]] + 1 ; r < B + 1 ; r++)
    {
        if (    f[k-1][r] + (y_til -x[itens[k-1]] ) < g[k] )
        {
            g[k] = f[k-1][r] + (y_til -x[itens[k-1]] );
                if (g[k] < g_best)
            {
                k_best = k;
                g_best = g[k];
            }   
        }
    }
    if (g[k] < y_til - 0.000001)
    {
        itens_adicioandos.push_back(itens[k-1]);
        soma_demandas_cobertura_auxiliar = soma_demandas_cobertura_auxiliar + d[itens[k - 1]];
                double demanda_minima = B + 1 - d[itens[k-1]];
        double soma_x = g[k] - (y_til -x[itens[k-1]] ); 
                              if (  (abs(soma_x) < 0.000001) && (demanda_minima < 0) ) 
            {
           itens_adicioandos.reverse();
        for (std::list<int>::iterator k1 = itens_adicioandos.begin(); k1 != itens_adicioandos.end(); k1++ ){
        vetor_cobertura_auxiliar.push_back(*k1);
        }
            corte_encontrado = true;
            cobertura_encontrada = true;
            return;
            }
        int x_pivo = k-1;
        int y_pivo;
        for (int r = demanda_minima; r < B+1; r++)
        {
           if ( abs(f[x_pivo][r] - soma_x) < 0.000001  )
           {
               y_pivo = r;
               break;
           } 
        }
        while (x_pivo >=1 )
        {
         if (x_pivo > 0)
         {
         if (  (abs (f[x_pivo - 1][y_pivo] - f[x_pivo][y_pivo])<0.000001   )  )
         {
            x_pivo = x_pivo - 1;
         }
         else
         {
               itens_adicioandos.push_back(itens[x_pivo-1]);
              soma_demandas_cobertura_auxiliar = soma_demandas_cobertura_auxiliar + d[itens[x_pivo - 1]];
              demanda_minima = demanda_minima - d[itens[x_pivo-1]];
              soma_x = soma_x - ( y_til - x[itens[x_pivo-1]]  ) ; 
                                            if (  (abs(soma_x)  < 0.000001) && (demanda_minima < 0) ) 
            {
            corte_encontrado = true;
            cobertura_encontrada = true;
                           itens_adicioandos.reverse();
        for (std::list<int>::iterator k1 = itens_adicioandos.begin(); k1 != itens_adicioandos.end(); k1++ ){
        vetor_cobertura_auxiliar.push_back(*k1);
        }
            return;
            }
              x_pivo = x_pivo - 1;
                for (int r = demanda_minima; r < B+1; r++)
                {
                if ( abs(f[x_pivo][r] - soma_x) < 0.000001  ) 
                {
                    y_pivo = r;
                break;
                } 
                }
         }
         }
        }// fim while
           itens_adicioandos.reverse();
        for (std::list<int>::iterator k1 = itens_adicioandos.begin(); k1 != itens_adicioandos.end(); k1++ ){
        vetor_cobertura_auxiliar.push_back(*k1);
        }
        corte_encontrado = true;
        cobertura_encontrada = true;
        return;   
    } // termina if o corte é encontrado
}   // fim para k = 1 ....
if (g_best < inf)
{
    cobertura_encontrada = true;
}
if ( (corte_encontrado == false) && (cobertura_encontrada == true)   )
{
   itens_adicioandos.push_back(itens[k_best-1]);
        soma_demandas_cobertura_auxiliar = soma_demandas_cobertura_auxiliar + d[itens[k_best - 1]];
                double demanda_minima = B + 1 - d[itens[k_best-1]];
        double soma_x = g_best- (y_til -x[itens[k_best-1]] ); 
           if (  (abs(soma_x) < 0.000001) && (demanda_minima < 0) ) 
            {
           itens_adicioandos.reverse();
        for (std::list<int>::iterator k1 = itens_adicioandos.begin(); k1 != itens_adicioandos.end(); k1++ ){
        vetor_cobertura_auxiliar.push_back(*k1);
        }
            return;
            }
        int x_pivo = k_best-1;
        int y_pivo;
        for (int r = demanda_minima; r < B+1; r++)
        {
           if ( abs(f[x_pivo][r] - soma_x) < 0.000001  )
           {
               y_pivo = r;
               break;
           } 
        }
        while (x_pivo >=1 )
        {
         if (x_pivo > 0)
         {
         if (  (abs (f[x_pivo - 1][y_pivo] - f[x_pivo][y_pivo])<0.000001   )  )
         {
            x_pivo = x_pivo - 1;
         }
         else
         {
               itens_adicioandos.push_back(itens[x_pivo-1]);
              soma_demandas_cobertura_auxiliar = soma_demandas_cobertura_auxiliar + d[itens[x_pivo - 1]];
              demanda_minima = demanda_minima - d[itens[x_pivo-1]];
              soma_x = soma_x - ( y_til - x[itens[x_pivo-1]]  ) ; 
                                            if (  (abs(soma_x)  < 0.000001) && (demanda_minima < 0) ) 
            {
                           itens_adicioandos.reverse();
        for (std::list<int>::iterator k1 = itens_adicioandos.begin(); k1 != itens_adicioandos.end(); k1++ ){
        vetor_cobertura_auxiliar.push_back(*k1);
        }
            return;
            }
              x_pivo = x_pivo - 1;
                for (int r = demanda_minima; r < B+1; r++)
                {
                if ( abs(f[x_pivo][r] - soma_x) < 0.000001  ) 
                {
                    y_pivo = r;
                break;
                } 
                }
         }
         }
        }// fim while
           itens_adicioandos.reverse();
        for (std::list<int>::iterator k1 = itens_adicioandos.begin(); k1 != itens_adicioandos.end(); k1++ ){
        vetor_cobertura_auxiliar.push_back(*k1);
        }
        return;
}
}  // fim else y_til = 1
} // fim void separacao exata programacao dinamica
void SolveKnapsack(double capacidade, vector<double> coef, vector<int> clientes_corte , vector<double> Pesos, double &valor_otimo, vector<int> &itens_otimo)
{

int n = clientes_corte.size();
int soma_weights = 0;
int soma_profits = 0;
int p[n];
for (int i = 0; i < n; i++)
{
    p[i] = coef[clientes_corte[i]];
    soma_profits = soma_profits + p[i];
}
int w[n];
for (int i = 0; i < n; i++)
{
    w[i] = Pesos[clientes_corte[i]];
    soma_weights = soma_weights + w[i];
}
int c = capacidade;
int x[n];
if (soma_weights <= c)
{
    valor_otimo = soma_profits;
    for (int i = 0; i < n; i++)
    {
       itens_otimo.push_back(i);
    }
}
else
{
valor_otimo = minknap(n, p, w, x ,c);
for (int i = 0; i < n; i++)
{
   if (x[i] > 0.5 ) itens_otimo.push_back(i);
}
}
}
void LCIsuperaditivefunction(vector<double> Pesos , vector<int> Itens, double B, vector<int> Cobertura, vector<double> &coef, vector<int> &clientes_corte, double &limitante, vector<double> &x, bool &corte_violado, double &y_til)
{
    funcao_ordena_decrescente(Cobertura,Pesos);
    vector<int> c_menos;
    vector<int> vetor_cobertura_auxiliar;
    limitante = Cobertura.size() - 1;
    double Soma = 0;
    vetor_cobertura_auxiliar.clear();
    vetor_cobertura_auxiliar.push_back(0); 
    for (int i = 0; i < Cobertura.size(); i++)
    {
      vetor_cobertura_auxiliar.push_back(Cobertura[i]);
      Soma = Soma + Pesos[Cobertura[i]];
    }
    int aux = Cobertura.size() ;
    double a_barra = Pesos[Cobertura[0]];
    double sigma = Soma - B;
 for (int k = 1; k < vetor_cobertura_auxiliar.size() - 1 ; k++)
    {
    double delta = a_barra - Pesos[vetor_cobertura_auxiliar[k+1]];
    if (k*delta < sigma)
    {
        a_barra = Pesos[vetor_cobertura_auxiliar[k+1]];
        sigma = sigma - k*delta;
    }
    else
    {
        a_barra = (double) a_barra - sigma/k;
        sigma = 0;
        break;
    }
}
if (sigma > 0)
{
    a_barra =  B/Cobertura.size() ;
}
list<int> lista_c_menos;
vector<int> c_mais;
for (int i = 0; i < Cobertura.size(); i++)
{
    if (Pesos[Cobertura[i]] <= a_barra + 0.0000000000001)
    {
        clientes_corte.push_back(Cobertura[i]);
        coef.push_back(1);
        lista_c_menos.push_back(Cobertura[i]);
    }
    else
    {
        c_mais.push_back(Cobertura[i]);
    }
}
vector<double> a_menos(Cobertura.size());
for (int i = 0; i < a_menos.size(); i++)
{
    if (Pesos[Cobertura[i]] <= a_barra)
    {
        a_menos[i] = Pesos[Cobertura[i]]; 
    }
    else
    {
        a_menos[i] = a_barra;
    }
    
}
int aux2 = Cobertura.size() + 1;
vector<double> S_menos(aux2);
S_menos[0] = 0;
for (int i = 1; i < aux2 ; i++)
{
    S_menos[i] = S_menos[i-1] + a_menos[i-1]; 
}
for (int i = 0; i < aux2; i++)
{
    if ( abs(S_menos[i] - round(S_menos[i]) ) < 0.00000001 ) S_menos[i] = round(S_menos[i]);
}
bool demanda_igual_a_barra = false;
for (int i = 0; i < Cobertura.size(); i++)
{
    if (abs(Pesos[Cobertura[i]] - a_barra) < 0.00001 )
    {
        demanda_igual_a_barra = true;
        break;
    }
}
list<double> lista_demandas_fortalecidas;
demanda_igual_a_barra = true;
if (demanda_igual_a_barra == true)
{
for (int h = 1; h < c_mais.size() ; h++)
{
    lista_demandas_fortalecidas.push_back(h*a_barra);
}
}
else
{
int inf = 0;
if (  (c_mais.size() % 2) == 0 )
{
 inf = ceil((double) c_mais.size()/2) + 1;
}
else
{
    inf = ceil((double) c_mais.size()/2);
}
for (int h = inf; h < c_mais.size() ; h++)
{
    lista_demandas_fortalecidas.push_back(h*a_barra);
}
}
vector<int> N1;
vector<int> N2;
for (int i = 0; i < Itens.size(); i++)
{
    if ( !esta_na_lista(lista_c_menos, Itens[i]) )
    {
     if (x[Itens[i]] > 0  )
        {
        N1.push_back(Itens[i]);
        }
        else
        {
        N2.push_back(Itens[i]);
        }
    }
}
if (demanda_igual_a_barra == true)
{
for (int i = 0; i < N1.size(); i++)
{
 if (Pesos[N1[i]] > B  )
 {
                coef.push_back(1);
            clientes_corte.push_back(N1[i]);
            continue;
 }
    int lambda = 0;
    bool pare2 = false;
    while (pare2 == false)
    {
         if (  ((Pesos[N1[i]] > S_menos[lambda]) && (Pesos[N1[i]] <= S_menos[lambda + 1])) )    
        {
            if (esta_na_lista_double(lista_demandas_fortalecidas, Pesos[N1[i]]))
            {
                lambda = lambda + 0.5;
            }
            if (lambda > 0)
            {
            coef.push_back(lambda);
            clientes_corte.push_back(N1[i]);
            }
            pare2 = true;
        }
        lambda = lambda + 1;
    }
}
}
else
{
for (int i = 0; i < N1.size() ; i++)
{
        if (Pesos[N1[i]] > B )
    {
            coef.push_back(1);
            clientes_corte.push_back(N1[i]);
            continue;
    }
    double lambda = 0;
    int lambda_anterior = 0;
    bool pare2 = false;
    while (pare2 == false)
    {
        if (  ((Pesos[N1[i]] > S_menos[lambda]) && (Pesos[N1[i]] <= S_menos[lambda + 1])) && !((Pesos[N1[i]] > S_menos[lambda+1]) && (Pesos[N1[i]] <= S_menos[lambda+2]))    )
        {
            if (esta_na_lista_double(lista_demandas_fortalecidas, Pesos[N1[i]]))
            {
                lambda = lambda + 1;
            }
            else
            {
                if (  ((c_mais.size() % 2) == 0 ) && (abs((c_mais.size()*a_barra)/2 - Pesos[N1[i]] ) < 0.00001 ) )   
                {
                    lambda = lambda + 0.5;
                }
            }
            if (lambda > 0)
            {
            coef.push_back(lambda);
            clientes_corte.push_back(N1[i]);
            }
            pare2 = true;
        }
        lambda = lambda + 1;
    }
}
}
if (  (y_til > 0.999999) && (verifica_corte_e_violado(Pesos , B, coef,clientes_corte, limitante,  x, y_til)  ) ) 
{
corte_violado = true;
if (demanda_igual_a_barra == true)
{
for (int i = 0; i < N2.size(); i++)
{
        if (Pesos[N2[i]] > B )
    {
            coef.push_back(1);
            clientes_corte.push_back(N2[i]);
            continue;
    }
    int lambda = 0;
    bool pare2 = false;
    while (pare2 == false)
    {
         if (  ((Pesos[N2[i]] > S_menos[lambda]) && (Pesos[N2[i]] <= S_menos[lambda + 1] + 0.000000001)) )    
        {
            if (esta_na_lista_double(lista_demandas_fortalecidas, Pesos[N2[i]]))
            {
                lambda = lambda + 0.5;
            }
            if (lambda >= 0)
            {
            coef.push_back(lambda);
            clientes_corte.push_back(N2[i]);
            }
            pare2 = true;
        }
        lambda = lambda + 1;
    }
}
}
else
{
for (int i = 0; i < N2.size() ; i++)
{
           if (Pesos[N2[i]] > B )
    {
            coef.push_back(1);
            clientes_corte.push_back(N2[i]);
            continue;
    }
    double lambda = 0;
    int lambda_anterior = 0;
    bool pare2 = false;
    while (pare2 == false)
    {
        if (  ((Pesos[N2[i]] > S_menos[lambda]) && (Pesos[N2[i]] <= S_menos[lambda + 1])) && !((Pesos[N2[i]] > S_menos[lambda+1]) && (Pesos[N2[i]] <= S_menos[lambda+2]))    )
        {
            if (esta_na_lista_double(lista_demandas_fortalecidas, Pesos[N2[i]]))
            {
                lambda = lambda + 1;
            }
            else
             {
                if (  ((c_mais.size() % 2) == 0 ) && (abs((c_mais.size()*a_barra)/2 - Pesos[N2[i]] ) < 0.00001 ) )   
                {
                    lambda = lambda + 0.5;
                }
            }
            if (lambda >= 0)
            {
            coef.push_back(lambda);
            clientes_corte.push_back(N2[i]);
            }
            pare2 = true;
        }
        lambda = lambda + 1;
    }
}
}
} // fim verifica se o corte é valido.
if (    y_til < 0.999999  ) 
{
if (demanda_igual_a_barra == true)
{
for (int i = 0; i < N2.size(); i++)
{
        if (Pesos[N2[i]] > B )
    {
            coef.push_back(1);
            clientes_corte.push_back(N2[i]);
            continue;
    }
    int lambda = 0;
    bool pare2 = false;
    while (pare2 == false)
    {
         if (  ((Pesos[N2[i]] > S_menos[lambda]) && (Pesos[N2[i]] <= S_menos[lambda + 1])) )    
        {
            if (esta_na_lista_double(lista_demandas_fortalecidas, Pesos[N2[i]]))
            {
                lambda = lambda + 0.5;
            }
            if (lambda > 0)
            {
            coef.push_back(lambda);
            clientes_corte.push_back(N2[i]);
            }
            pare2 = true;
        }
        lambda = lambda + 1;
    }
}
}
else
{
for (int i = 0; i < N2.size() ; i++)
{
           if (Pesos[N2[i]] > B )
    {
            coef.push_back(1);
            clientes_corte.push_back(N2[i]);
            continue;
    }
    double lambda = 0;
    int lambda_anterior = 0;
    bool pare2 = false;
    while (pare2 == false)
    {
        if (  ((Pesos[N2[i]] > S_menos[lambda]) && (Pesos[N2[i]] <= S_menos[lambda + 1])) && !((Pesos[N2[i]] > S_menos[lambda+1]) && (Pesos[N2[i]] <= S_menos[lambda+2]))    )
        {
            if (esta_na_lista_double(lista_demandas_fortalecidas, Pesos[N2[i]]))
            {
                lambda = lambda + 1;
            }
            else
             {
                if (  ((c_mais.size() % 2) == 0 ) && (abs((c_mais.size()*a_barra)/2 - Pesos[N2[i]] ) < 0.00001 ) )   
                {
                    lambda = lambda + 0.5;
                }
            }
            if (lambda > 0)
            {
            coef.push_back(lambda);
            clientes_corte.push_back(N2[i]);
            }
            pare2 = true;
        }
        lambda = lambda + 1;
    }
}
}

if (verifica_corte_e_violado(Pesos , B, coef,clientes_corte, limitante,  x, y_til) )
{
    corte_violado = true;
}
} // fim verifica se o corte é valido.
} // fim void

bool verifica_corte_e_valido_usando_solver(vector<double> Pesos , double B, vector<double> coef, vector<int> clientes_corte, double limitante)
{
vector<int> itens_mochila;
vector<double> pesos_mochila;
vector<double> valores_mochila;
int capacidade = B;

for (int i = 0; i < clientes_corte.size(); i++)
{
    itens_mochila.push_back(i);
    pesos_mochila.push_back(Pesos[clientes_corte[i]]);
    valores_mochila.push_back(coef[i]);
}

IloEnv env_mochila;
IloModel mod_mochila(env_mochila);
IloCplex cplex_mochila(mod_mochila);
IloNumVarArray x_mochila(env_mochila, itens_mochila.size(), 0, 1, ILOBOOL); 
// objective *****************************************************
IloExpr expfo_mochila(env_mochila);
for (int i = 0; i <  itens_mochila.size(); i++){
	expfo_mochila += coef[i] * x_mochila[i];

}  
IloAdd(mod_mochila, IloMaximize(env_mochila, expfo_mochila));
expfo_mochila.end();

IloExpr rest_mochila(env_mochila);
for (int i = 0; i < itens_mochila.size(); i++){
			rest_mochila += x_mochila[i]*pesos_mochila[i];
}
mod_mochila.add( rest_mochila <= B)   ;
cplex_mochila.setWarning(env_mochila.getNullStream()); // Eliminar warnings
cplex_mochila.setOut(env_mochila.getNullStream()); /// Eliminar os logs do solver
cplex_mochila.setParam(IloCplex::Threads,1);
cplex_mochila.setParam(IloCplex::Param::Parallel,1);
cplex_mochila.solve();
double mochila = cplex_mochila.getObjValue();  
if (limitante < mochila )
{
    return false;
}
else
{
    return true;
}
} //fim

bool verifica_corte_e_violado(vector<double> Pesos , double B, vector<double> coef, vector<int> clientes_corte, double limitante, vector<double>  x, double &y_til)
{
double soma = 0;
for (int i = 0; i < clientes_corte.size(); i++)
{
    soma = soma + x[clientes_corte[i]]*coef[i];
}
if (soma > limitante*y_til)
{
    return true;
}
else
{
    return false;
}
}

bool verifica_corte_e_violadoB(vector<double> Pesos , double B, vector<double> coef, vector<int> clientes_corte, double limitante, vector<double>  x)
{
double soma = 0;
for (int i = 0; i < clientes_corte.size(); i++)
{
    soma = soma + x[clientes_corte[i]]*coef[i];
}
if (soma > limitante)
{
    return true;
}
else
{
    return false;
}





}

void LCIsequentiallifting(vector<double> Pesos , vector<int> Itens, double B, vector<int> Cobertura, vector<double> &coef, vector<int> &clientes_corte, double &limitante, vector<double> x, double &y_til)
{
  list<int> lista_cobertura;
  list<int> lista_D;
  for (int i = 0; i < Cobertura.size(); i++)
  {
      lista_cobertura.push_back(Cobertura[i]);
  }
vector<int> CmenosD;
double soma_das_demandasD = 0;
   for (int i = 0; i < Cobertura.size(); i++)
    {
      if (x[Cobertura[i]] == 1 )
      {
         lista_D.push_back(Cobertura[i]);
         
         soma_das_demandasD = soma_das_demandasD + Pesos[Cobertura[i]];
      }
      else
      {
          CmenosD.push_back(Cobertura[i]);
      }
    }
  vector<int> N1;  //posições que não estão na cobertura mas o valor de x é positivo
  vector<int> N2; //posições que não estão na cobertura mas o valor de x é zero
  int maior_demanda_em_N1 = 0;
  for (int i = 0; i < Itens.size(); i++)
  {
     if (!esta_na_lista (lista_cobertura, Itens[i]) ) 
     {
      if ((x[Itens[i]] > 0  ) )
       {
            N1.push_back(i);  
            if (Pesos[i] > maior_demanda_em_N1)
            {
                maior_demanda_em_N1 = Pesos[i];
            } 
       }
       else
       {
            N2.push_back(i);
       }
     }
  }

double capacidade = B - soma_das_demandasD;
while (maior_demanda_em_N1 > capacidade)
{
  capacidade = capacidade + Pesos[lista_D.front()];
  CmenosD.push_back(lista_D.front());
  if (Pesos[lista_D.front()] > maior_demanda_em_N1)
  {
      maior_demanda_em_N1 = Pesos[lista_D.front()];
  }  
}
vector<int> D = cria_vetor_copia_de_lista(lista_D);
limitante = CmenosD.size() - 1;
clientes_corte = CmenosD;
for (int i = 0; i < CmenosD.size(); i++)
{
    coef.push_back(1);
}
clock_t tempo1 = clock();
if (N1.size() > 0)
{
for (int i = 0; i < N1.size(); i++)
{
vector<int> itens_mochila;
vector<double> pesos_mochila;
vector<double> valores_mochila;
for (int i = 0; i < clientes_corte.size(); i++)
{
    itens_mochila.push_back(i);
    pesos_mochila.push_back(Pesos[clientes_corte[i]]);
    valores_mochila.push_back(coef[i]);
}
capacidade = capacidade -  Pesos[N1[i]];
double otimo_mochila;
vector<int> itens_otimo;
SolveKnapsack( capacidade, valores_mochila, itens_mochila , pesos_mochila, otimo_mochila, itens_otimo) ;
int mochila = otimo_mochila;
int coef_aux = limitante - mochila;
if (coef_aux > 0)
{
clientes_corte.push_back(N1[i]);
coef.push_back(coef_aux); 
}
capacidade = capacidade +  Pesos[N1[i]];
}
}
clock_t fim_tempo1 = clock();
clock_t tempo2 = clock();
if (D.size() > 0)
{
for (int i = 0; i < D.size(); i++)
{
vector<int> itens_mochila;
vector<double> pesos_mochila;
vector<double> valores_mochila;
for (int i2 = 0; i2 < clientes_corte.size(); i2++)
{
    itens_mochila.push_back(i2);
    pesos_mochila.push_back(Pesos[clientes_corte[i2]]);
    valores_mochila.push_back(coef[i2]);
}
capacidade = capacidade +  Pesos[D[i]];
double otimo_mochila;
vector<int> itens_otimo;
SolveKnapsack( capacidade, valores_mochila, itens_mochila , pesos_mochila, otimo_mochila, itens_otimo) ;
int mochila = otimo_mochila;
clientes_corte.push_back(D[i]);
int coef_aux = mochila - limitante;
coef.push_back(coef_aux); 
limitante = limitante + coef_aux;
}
}
clock_t fim_tempo2 = clock();
clock_t tempo3 = clock();
if ( ( 1> 0) && (N2.size() > 0) &&  ( verifica_corte_e_violado(Pesos ,B, coef, clientes_corte,limitante, x, y_til) )   )
{
bool opcao1 = false;  //estrategia ECI
bool opcao2 = false; // lifting sequencial
if (opcao1 == true)
{
  int maior_peso = 0;
  int indice_cliente_de_maior_peso;
for (int i = 0; i < clientes_corte.size(); i++)
{
    if (Pesos[clientes_corte[i]] > maior_peso )
    {
        maior_peso = Pesos[clientes_corte[i]];
        indice_cliente_de_maior_peso = i; 
    }
}
for (int i = 0; i < N2.size(); i++)
{
    if (Pesos[N2[i]]>= maior_peso)
    {
        clientes_corte.push_back(N2[i]);
        coef.push_back(coef[indice_cliente_de_maior_peso]);
    }
}  
}
else
{
for (int i = 0; i < N2.size(); i++)
{
vector<int> itens_mochila;
vector<double> pesos_mochila;
vector<double> valores_mochila;
for (int i2 = 0; i2 < clientes_corte.size(); i2++)
{
    itens_mochila.push_back(i2);
    pesos_mochila.push_back(Pesos[clientes_corte[i2]]);
    valores_mochila.push_back(coef[i2]);
}
capacidade = capacidade - Pesos[N2[i]]; // qual é o valor máximo considerando as variáveis do corte caso a variável presente entre no corte com valor 1.  
if(capacidade < 0) continue;
vector<int> itens_otimo;
double otimo_mochila = 0;
SolveKnapsack( capacidade, valores_mochila, itens_mochila , pesos_mochila, otimo_mochila, itens_otimo) ;
int coef_aux = round(limitante - otimo_mochila);
if (coef_aux != 0)
{
coef.push_back(coef_aux); 
clientes_corte.push_back(N2[i]);
}
capacidade = capacidade +  Pesos[N2[i]];   //volta a capacidade que estava antes
}
}

} //fim se N2.size() > 0
clock_t fim_tempo3 = clock();
}   //fim constroi lci geral

void LPIsuperaditivefunction(vector<double> Pesos , vector<int> Itens, double B, vector<int> vetor_pack, list<int> pack, vector<double> &alfa, vector<double> &lambda, vector<int> &vetor_clientes_restantes, double &limitante, vector<double> x, bool &corte_violado, double soma_demandas_pack)
{
    //vector<int> vetor_clientes_restantes;
    double residual = B - soma_demandas_pack;
    for (int i = 0; i < Itens.size(); i++)
    {
    if (!esta_na_lista(pack, Itens[i]) && (Pesos[Itens[i]]) > residual )
    {
        vetor_clientes_restantes.push_back(Itens[i]);
    }
    }
    funcao_ordena_decrescente(vetor_clientes_restantes,Pesos);
    for (int i = 0; i < vetor_clientes_restantes.size(); i++)
    {
    int aux = Pesos[vetor_clientes_restantes[i]] - residual;
    if (aux >= 0)
    {
        alfa.push_back(aux);
    }
    else
    {
        alfa.push_back(0);
    }
    }
    vector<int> A(vetor_clientes_restantes.size()+ 1);
    A[0] = 0;
    for (int i = 1; i < vetor_clientes_restantes.size() + 1; i++)
    {
    A[i] = A[i-1] - Pesos[vetor_clientes_restantes[i-1]];
    }
    for (int k = 0; k < vetor_pack.size(); k++)
    {
    int i = 0;
    bool encontrou = false;
    while (encontrou == false)
    {
       if (  (-Pesos[vetor_pack[k]] >= (A[i+1] + residual))  && (-Pesos[vetor_pack[k]] <= A[i] )  )
       {
            lambda.push_back(i*residual - Pesos[vetor_pack[k]]);
            encontrou = true;
            break;
       }
       if ( (-Pesos[vetor_pack[k]] >= A[i]) && (-Pesos[vetor_pack[k]] <= (A[i] + residual)) )
       {
            lambda.push_back(i*residual + A[i]);
            encontrou = true;
            break;
       }
              if ( -Pesos[vetor_pack[k]] <= (A[vetor_clientes_restantes.size()] + residual)   )
       {
            lambda.push_back(vetor_clientes_restantes.size()*residual + A[vetor_clientes_restantes.size()] );
            encontrou = true;
            break;
       }
    i = i + 1;    
    }
    
}

}

bool inteiroEstaNoVetor( vector<int> &vetor, int valor) {
    for (int elemento : vetor) {
        if (elemento == valor) {
            return true;
        }
    }
    return false;
}

bool vetorIdenticoNaMatriz( vector<vector<int>> matriz, vector<int> vetor, int k) {
    
    for (int i = 0; i <= k; ++i) {
        if (matriz[i] == vetor) {
            //cout << "é igual a linha:" << i << endl;
            return true;
        }
    }
    return false;
    //cout << endl;
}

int indice_de_valor_minimo(const std::vector<int>& vetor1, const std::vector<int>& vetor2) {
    if (vetor1.empty() || vetor2.empty()) {
        return -1;  // Retorna -1 se os vetores estiverem vazios
    }

    int minIndex = 0;
    int minValue = vetor2[0];

    for (int i = 1; i < vetor1.size(); i++) {
        if (vetor2[i] < minValue) {
            minValue = vetor2[i];
            minIndex = i;
        }
    }

    return minIndex;
}

int valor_minimo(const std::vector<int>& vetor1, const std::vector<int>& vetor2) {
    if (vetor1.empty() || vetor2.empty()) {
        return -1;  // Retorna -1 se os vetores estiverem vazios
    }

    int minIndex = 0;
    int minValue = vetor2[0];

    for (int i = 1; i < vetor1.size(); i++) {
        if (vetor2[i] < minValue) {
            minValue = vetor2[i];
            minIndex = i;
        }
    }

    return minValue;
}

void imprime_vetor(vector<int> vetor) 
{
    for (int i = 0; i < vetor.size(); i++)
    {
        cout << vetor[i] << " ";
    }
    
}

vector<int> cria_vetor_copia_de_lista(const std::list<int>& myList) {
    std::vector<int> result(myList.begin(), myList.end());
    return result;
}

bool esta_na_lista(const std::list<int>& lista, int numero) {    //verifica se um número está em uma lista
    for (const int& elemento : lista) {
        if (elemento == numero) {
            return true;
        }
    }
    return false;
}

bool esta_na_lista_double(const std::list<double>& lista, int numero) {    //verifica se um número está em uma lista
    for (const double& elemento : lista) {
        //cout << "numero:" << numero << endl;
        //cout << "elemento da lista:" << elemento << endl;
        //cout << abs(elemento - numero) << endl; 
        if ( abs(elemento - numero) < 0.00000001 ) {
            
            return true;
        }
    }
    return false;
}

void cria_copia_lista (list<int> &lista_nova, list<int> lista)
{
lista_nova.clear();
for (std::list<int>::iterator k = lista.begin(); k != lista.end(); k++ ){
lista_nova.push_back(*k);
}
}

void cria_copia_vetor_int(vector<int> &vetor_novo, vector<int> vetor)
{
for (int i = 0; i < vetor.size(); i++)
{
   vetor_novo[i] = vetor[i];
}

}

void confere_solucao_s(int m, int n, vector<vector<double>> c,  vector<double> f, vector<list<int>>s, vector<double> d, list<int> abertas, list<int> fechadas, vector<double> pen, vector<double> p)
{
double fo1 = 0; //custos abertura facilidades
double fo2 = 0; // custos atendimento clientes
int soma_clientes = 0;   
int conta = 0; 
vector<double> penalidade(m);
for (int j = 0; j < m; j++)
{
    penalidade[j] = -p[j];
}

for (int j = 0; j < m; j++)
{
     
   s[j].sort();
   s[j].unique();
   if (!s[j].empty())
   {
        fo1 = fo1 + f[j]; 
        for (std::list<int>::iterator k1 = s[j].begin(); k1 != s[j].end(); k1++ ){
            soma_clientes = soma_clientes +1 ;
            fo2 = fo2 + c[*k1][j]*d[*k1];
            penalidade[j] = penalidade[j] + d[*k1];
            }
   }
   if (penalidade[j] > 0)
   {
    conta = conta + 1;
   }
}

double fo3;
  for (std::list<int>::iterator k = abertas.begin(); k != abertas.end(); k++ ){
  fo3 = fo3 + f[*k];
  }
int erro = 0;
for (int j = 0; j < m; j++)
{
    if (pen[j] != penalidade[j])
    {
        erro = 1;
        break;
    }
}



cout << fo1 + fo2 << " " << soma_clientes <<  " " << fo1 << " " << fo2 << " " << fechadas.size() + abertas.size() << " " << erro << endl;



}

void confere_solucao_atende (int m, int n, vector<vector<double>> c,  vector<double> f, vector<int> atende, vector<double> d, vector<double> pen, double &fo_confere)
{

double fo1_confere = 0 ;
double fo2_confere = 0;
vector<int> abertas_auxiliar(m);
for (int i = 0; i < abertas_auxiliar.size(); i++)
{
    abertas_auxiliar[i] = 0;
}


for (int i = 0; i < n; i++)
{
    for (int j = 0; j < m; j++)
    {
        if (atende[i] == j)
        {
            fo1_confere = fo1_confere + c[i][j]*d[i];
            abertas_auxiliar[j] = 1;
        }
    }
    
}

for (int j = 0; j < m; j++)
{
    fo2_confere = fo2_confere + abertas_auxiliar[j]*f[j] + pow(10,5)*max(pen[j],0.0);
}


/*
cout << "abertas atende:" << endl;
for (int j = 0; j < m; j++)
{
    if (abertas_auxiliar[j] == 1)
    {
        cout << j << " ";
    }
}
*/
cout << endl;
cout << setprecision(16) << "soma dos custos fixos de abertura: " << fo2_confere << endl;
cout << setprecision(16)  << "soma dos custos atendimento: " << fo1_confere << endl;
fo_confere = fo1_confere + fo2_confere;
cout << setprecision(16) << " fo confere: " << fo1_confere + fo2_confere << " ";
}

bool verifica_factivel( vector<double> pen)
{
    for (int j = 0 ; j < pen.size() ; j++)
    {
        if (pen[j] > 0)
        {
        return false;
        }
    }
    return true;
    
}

void funcao_ordena_decrescente(vector<int>&vetor_clientes, vector<double> d)
{
//cout << ".";


    auto compararPorD = [&d](int a, int b) {
        return d[a] > d[b]; // Ordenar em ordem decrescente com base em d
    };
    std::sort(vetor_clientes.begin(), vetor_clientes.end(), compararPorD);

}

void funcao_ordena_crescente(vector<int>&vetor_clientes, vector<double> d)
{



    auto compararPorD = [&d](int a, int b) {
        return d[a] < d[b]; // Ordenar em ordem decrescente com base em d
    };
    std::sort(vetor_clientes.begin(), vetor_clientes.end(), compararPorD);

}

void funcao_ordena_decrescente2(vector<int>& vetor_clientes, const vector<double>& x, const vector<double>& d) {

    // Função lambda para comparação
    auto comparar = [&x, &d](int a, int b) {
        if (x[a] > x[b]) {
            return true; // a deve vir antes de b
        } else if (x[a] < x[b]) {
            return false; // b deve vir antes de a
        } else {
            // x[a] == x[b], verificar por d em ordem decrescente
            return d[a] > d[b];
        }
    };

    // Usando a função lambda no sort
    sort(vetor_clientes.begin(), vetor_clientes.end(), comparar);
}

void ordena_clientes_crescente_demanda(vector<int>&vetor_clientes, vector<double> d)
{
 for (int i = 0; i < vetor_clientes.size()- 1; i++)
    {
    int menor_demanda = i;
    for (int j = i; j < vetor_clientes.size(); j++)
    {
       if (d[vetor_clientes[j]] < d[vetor_clientes[menor_demanda]])
       {
           menor_demanda = j;
       }
    }
    if (menor_demanda != i)
    {
        int temp = vetor_clientes[i];
        vetor_clientes[i] = vetor_clientes[menor_demanda];
        vetor_clientes[menor_demanda] = temp;
    }
    
    }

}

void imprime_vetor_double(vector<double> vetor)
{
    for (int i = 0; i < vetor.size(); i++)
    {
       cout << setprecision(12) << vetor[i] << " ";
    }
    cout << endl;
}

void imprime_lista(list<int> lista)
{
for (  std::list<int>::iterator k = lista.begin() ; k !=  lista.end() ; k++)
{
cout << *k << " ";
}

}

void imprime_lista_double(const std::list<double>& minhaLista) {
    for (const double& valor : minhaLista) {
        std::cout << valor << " ";
    }
    std::cout << std::endl;
}

bool verifica_factivel_restrito(vector<int> vetor_facilidades, vector<double> pen)
{
    for (int j = 0 ; j < vetor_facilidades.size() ; j++)
    {
        if (pen[vetor_facilidades[j]] > 0)
        {
        return false;
        }
    }
    return true;
    
}

void embaralhar( vector<int> &ordem, int numero_vizinhancas)
{
	

    for (int i = 0; i < numero_vizinhancas; i++)
{
    ordem[i] = i;
}
    
    for (int i = 0; i < numero_vizinhancas; i++)
	{
		int r = rand() % numero_vizinhancas;

		int temp = ordem[i];
		ordem[i] = ordem[r];
		ordem[r] = temp;
	}
}

int verifica_lista_tabu (vector<vector<int>> lista_tabu, vector<int> vetor_facilidades_proximas1, int Q, int tamanho_lista_tabu, int ultima_melhora)
{


   for (int j = ultima_melhora; j < tamanho_lista_tabu + 1; j++)   //para cada vetor da lista tabu
   {



   int conta = 0; //compara posições do cluster corrente com as posições do vetor a lista tabu


   for (int i1 = 0; i1< vetor_facilidades_proximas1.size(); i1++) // verifica para o cluster corrente o numero de posicoes diferentes
       {

         //verifica se vetor_facilidades_proximas 1 [i1] está
         bool encontrou = false;   //verifica se a posicao i do cluster está na lista tabu posicao j
         for (int i2 = 0; i2 < Q; i2++)   //percorre as posições do cluster que esta na lista tabu
         {
           if (vetor_facilidades_proximas1[i1] == lista_tabu[j][i2])
           {
               encontrou = true;  //encontrou pelo menos uma igual
               break;
           }
         }

        if (encontrou == false)
         {
             conta = conta + 1;  //quantidade de diferentes
         }

        }

        if (conta == 0)  // o cluster está na lista
        {

          return 1;
        }


       /*

        cout << "conta posicoes diferentes cluster posicao " << j << " conta = : " <<  conta << endl;
        if ( conta <= menor_diferenca_permitida)   // se o numero de posicoes diferentes é menor ou igual ao permitido para
        {
            return 0;
        }
        */
     }

    return 0;
} //fim void

void acrescenta_lista_tabu (vector<vector<int>> &lista_tabu, vector<int> cluster, int Q, int tamanho_lista_tabu)
{
    //cout << "acrescenta lista tabu" << endl;
    //cout << "tamanho cluster: " << cluster.size() << endl;
    //cout << "valor de Q : " << Q << endl;
    //cout << "tamanho lista tabu : " << tamanho_lista_tabu << endl;

for (int i = 0; i < Q; i++)
{
    lista_tabu[tamanho_lista_tabu][i] = cluster[i];
}
  //cout << "fim lista tabu" << endl;
}

int verifica_lista_tabu_facilidades (vector<vector<int>> lista_tabu_facilidades, vector<int> vetor_facilidades_proximas1, int Q, int tamanho_lista_tabu_facilidades, int ultima_melhora)
{


   for (int j = ultima_melhora; j < tamanho_lista_tabu_facilidades + 1; j++)   //para cada vetor da lista tabu
   {


       
   

   int conta = 0; //compara posições do cluster corrente com as posições do vetor a lista tabu


   for (int i1 = 0; i1< vetor_facilidades_proximas1.size(); i1++) // verifica para o cluster corrente o numero de posicoes diferentes
       {

         //verifica se vetor_facilidades_proximas 1 [i1] está
         bool encontrou = false;   //verifica se a posicao i do cluster está na lista tabu posicao j
         for (int i2 = 0; i2 < Q; i2++)   //percorre as posições do cluster que esta na lista tabu
         {
           if (vetor_facilidades_proximas1[i1] == lista_tabu_facilidades[j][i2])
           {
               encontrou = true;  //encontrou pelo menos uma igual
               break;
           }
         }

        if (encontrou == false)
         {
             conta = conta + 1;  //quantidade de diferentes
         }

        }

        if (conta == 0)  // o cluster está na lista
        {

          return 1;
        }


       /*

        cout << "conta posicoes diferentes cluster posicao " << j << " conta = : " <<  conta << endl;
        if ( conta <= menor_diferenca_permitida)   // se o numero de posicoes diferentes é menor ou igual ao permitido para
        {
            return 0;
        }
        */
     }

    return 0;
} //fim void

void acrescenta_lista_tabu_facilidades (vector<vector<int>> &lista_tabu_facilidades, vector<int> cluster, int Q, int tamanho_lista_tabu_facilidades)
{
    //cout << "acrescenta lista tabu" << endl;
    //cout << "tamanho cluster: " << cluster.size() << endl;
    //cout << "valor de Q : " << Q << endl;
    //cout << "tamanho lista tabu : " << tamanho_lista_tabu << endl;

for (int i = 0; i < cluster.size(); i++)
{
    lista_tabu_facilidades[tamanho_lista_tabu_facilidades].push_back(cluster[i]);
}
  //cout << "fim lista tabu" << endl;
}

int verifica_lista_tabu_clientes (vector<vector<int>> lista_tabu_clientes, vector<int> vetor_conjunto_clientes, int Q, int tamanho_lista_tabu_clientes, int ultima_melhora)
{


   for (int j = ultima_melhora; j < tamanho_lista_tabu_clientes + 1; j++)   //para cada vetor da lista tabu
   {


   if (vetor_conjunto_clientes.size() != lista_tabu_clientes[j].size())
   {
       continue;
   }    
   int conta = 0; //compara posições do cluster corrente com as posições do vetor a lista tabu


   for (int i1 = 0; i1< vetor_conjunto_clientes.size(); i1++) // verifica para o cluster corrente o numero de posicoes diferentes
       {

         
         //verifica se vetor_facilidades_proximas 1 [i1] está
         bool encontrou = false;   //verifica se a posicao i do cluster está na lista tabu posicao j
         for (int i2 = 0; i2 < lista_tabu_clientes[j].size(); i2++)   //percorre as posições do cluster que esta na lista tabu
         {
           if (vetor_conjunto_clientes[i1] == lista_tabu_clientes[j][i2])
           {
               encontrou = true;  //encontrou pelo menos uma igual
               break;
           }
         }

        if (encontrou == false)
         {
             conta = conta + 1;  //quantidade de diferentes
         }

        }

        if (conta == 0)  // o cluster está na lista
        {

          return 1;
        }


       /*

        cout << "conta posicoes diferentes cluster posicao " << j << " conta = : " <<  conta << endl;
        if ( conta <= menor_diferenca_permitida)   // se o numero de posicoes diferentes é menor ou igual ao permitido para
        {
            return 0;
        }
        */
     }

    return 0;
} //fim void

void acrescenta_lista_tabu_clientes (vector<vector<int>> &lista_tabu_clientes, vector<int> vetor_conjunto_clientes, int tamanho_lista_clientes)
{
    //cout << "acrescenta lista tabu" << endl;
    //cout << "tamanho cluster: " << cluster.size() << endl;
    //cout << "valor de Q : " << Q << endl;
    //cout << "tamanho lista tabu : " << tamanho_lista_tabu << endl;

for (int i = 0; i < vetor_conjunto_clientes.size(); i++)
{
    lista_tabu_clientes[tamanho_lista_clientes].push_back(vetor_conjunto_clientes[i]);
}
  //cout << "fim lista tabu" << endl;
}


int verifica_lista_tabu_novo(vector<vector<int>> lista_tabu, vector<int> vetor, int inicio_da_lista, int final_da_lista) 
{


   for (int j = inicio_da_lista; j < final_da_lista + 1; j++)   //para cada vetor da lista tabu
   {


    if (vetor.size() != lista_tabu[j].size())
    {
           continue;
    }    


   bool igual = true; 
   for (int i1 = 0; i1< vetor.size(); i1++) // verifica para o cluster corrente o numero de posicoes diferentes
       {
            if (vetor[i1] != lista_tabu[j][i1])
            {
                igual = false;
            }
       }
    if (igual == false) 
    {
    continue;
    }
    else
    {
        return 1;
    }
   }
}

void verifica_se_problema_esparso_e_gap(list<int> locais_promissores, bool &problema_esparso_e_gap, vector<int> vetor_conjunto_clientes, vector<double> p, double soma_demandas, vector<int> vetor_facilidades_proximas1, double &soma_capacidades, vector<double> pen)
{

//verifica se há alguma facilidade vazia;
for (int i = 0; i < vetor_facilidades_proximas1.size(); i++)
{
if (pen[vetor_facilidades_proximas1[i]] == -p[vetor_facilidades_proximas1[i]] )
{
    //cout << "há facilidade vazia" << endl;
    problema_esparso_e_gap = false;
    return;
}
}


int indice_aberta_de_menor_capacidade = 0;
for (int j = 0; j < vetor_facilidades_proximas1.size(); j++)
{
   if ( p[vetor_facilidades_proximas1[j]] < p[vetor_facilidades_proximas1[indice_aberta_de_menor_capacidade]]  )
   {
    indice_aberta_de_menor_capacidade = j;
   }
}

soma_capacidades = 0;
for (int j = 0; j < vetor_facilidades_proximas1.size(); j++)
{
    soma_capacidades = soma_capacidades + p[vetor_facilidades_proximas1[j]];
}

if((soma_capacidades - soma_demandas) < p[vetor_facilidades_proximas1[indice_aberta_de_menor_capacidade]]  )
{
    //out << "há facilidade vazia" << endl;
    problema_esparso_e_gap = true;
   // cout << "problema esparso é gap" << endl;
    return;
}


}




