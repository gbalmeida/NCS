#include <iostream>
#include <iomanip>
#include <fstream>
#include <algorithm>
#include <math.h>
#include <vector>
#include <list>
#include <iterator>
#include <stdlib.h>
#include <ctime>
#include <time.h>
#include <chrono>
#include <ilcplex/ilocplex.h>
#include "Functions.h"
#include <ctime>

using namespace std;

ILOSTLBEGIN 

int main (int argc, char *argv[]) {
	try {    

srand (time(NULL));
        clock_t inicio_CPU, fim_CPU;
        inicio_CPU = clock();

// cutting-plane torerance
double tolerancia_pc = 0.000001; 
int bigM = 10e6;
// time limit to cutting-plane procedure
int tempo_maximo_pc = 3600;   
// time limit to solve linear partial relaxation 
double tempo_limite_problema_estrategia_integralidade_parcial = 3600; 
// max iter to generate fenchel cuts
int max_iter_fc = 50;
// time limit to NCS
int tempo_limite_total = 28800;

////////// select cutting-planes procedures generation //////////////////////

bool aplica_vubs = true;
bool aplica_CIgreedyA = true;
bool aplica_CIgreedyB = true;
bool aplica_CIgreedyNew = true;
bool aplica_CIexact = true;
bool aplica_LCIsuperaditivefunction = true;
bool aplica_LCIsequentiallifting = true;
bool aplica_LPIsuperaditivefunction = false;
bool aplica_FCColumngeneration = true;

///////// Heuristic parameters //////// 
int Qmin = 6;
int R = 4;
int U = 10;
int tempo_limite_heuristica = 3600;
///////////////////////////////////////

//////////// data structures for storing generated cuts ////////
list<vector<double>> coeficientes_cortes_x;
list<vector<int>> clientes_cortes_x;
vector<double> limitantes_cortes_x;
vector<int> variavel_y_associada_ao_corte;
list<vector<double>> coeficientes_cortes_y;
list<vector<int>> facilidades_cortes_y;
vector<double> limitantes_cortes_y;
////////////////////////////////////////////////////////////////



int m, n;
vector<double> p, f, d, pen; 
vector<vector<double>> c;
double soma_demandas, capacidade_total;
bool coeficientes_inteiros = false;

/// read instances ///////////

    Read_instances(argv[1], m, n, p, f, d, c, pen, soma_demandas, capacidade_total, coeficientes_inteiros);

double demanda_total = soma_demandas;


////////////// auxiliary variables////////////
//////////////////////////////////////////////
vector<list<int>> s(m);
double fo;
double tempo_total = 0; 
double tempo_matheuristica = 0;
double tempo_resolve_problema_denso = 0;
double tempo_primeiro_UB = 0;
bool solucao_inteira = false;
list<int> abertas, fechadas, locais_promissores;
vector<int> atende(n);
double soma_capacidades_locais_promissores = 0;
double LB_modelo_otimo; 
double lb_esparso = pow(10,10);
double primeiro_lb;
double primeiro_UB = 0;
double tempo_primeiro_lb;
int iteracao_que_encontrou_otimo = 0;
int quantidade_iteracoes = 0;
double otimo;
int maxiter = 1000;
double LB = 0;
double UB = 10e10;
bool primeira_solucao_factivel_encontrada = false;
bool factivel = true;
double tempo_relaxacao_no_raiz;
double limitante_inferior_no_raiz;
bool roda_heuristica = true;
bool imprime_detalhes = true;
int tempo_limite_problema_esparso;
vector<int> melhor_atende(n);
int iter = 0;  
bool solucao_heuristica = false;  
/////////////////////////////////////////////
/////////////////////////////////////////////





////////////////// linear relaxation model ///////////////////

IloEnv env_denso;
IloModel mod_denso(env_denso);
IloCplex cplex_denso(mod_denso);
IloNumVarArray y_denso(env_denso, m, 0, 1, ILOFLOAT); 
IloArray<IloNumVarArray> x_denso(env_denso, n);
for(int i = 0; i < n; i++) {  
			 x_denso[i] = IloNumVarArray(env_denso, m, 0, 1, ILOFLOAT);
}
list<int> locais_nao_promissores; // conjunto de locais nao promissores
IloConstraintArray restricoes_fcs(env_denso); // conjunto de restricoes do tipo fenchel cuts
cplex_denso.setParam(IloCplex::ParallelMode, 1 );
   vector<int> Itens(n);
        for (int i = 0; i < n; i++)
        {
          Itens[i] = i; 
        }
           vector<int> vetor_facilidades_proximas1(m);
        for (int i = 0; i < m; i++)
        {
          vetor_facilidades_proximas1[i] = i; 
        }
vector<double> y_copia(m);
vector<vector<double>> x_copia(n, vector<double>(m) )   ;

////////////////// run cutting-plane procedure ///////////////////
//////////////////////////////////////////////////////////////////

CuttingPlane(m,n, f, p,c , LB ,soma_demandas, d, y_copia, x_copia, capacidade_total, demanda_total,  fim_CPU, inicio_CPU, quantidade_iteracoes, iteracao_que_encontrou_otimo, tempo_resolve_problema_denso, tolerancia_pc,  bigM, aplica_LCIsuperaditivefunction, aplica_LPIsuperaditivefunction, aplica_LCIsequentiallifting, aplica_FCColumngeneration,  coeficientes_cortes_x, clientes_cortes_x, limitantes_cortes_x, variavel_y_associada_ao_corte, coeficientes_cortes_y,facilidades_cortes_y, limitantes_cortes_y,  max_iter_fc, aplica_CIgreedyA, aplica_CIgreedyB, aplica_CIexact, aplica_vubs, aplica_CIgreedyNew, tempo_maximo_pc,  tempo_limite_problema_estrategia_integralidade_parcial, env_denso, mod_denso, cplex_denso, y_denso,  x_denso, solucao_inteira, restricoes_fcs);

if (solucao_inteira == true)
{
UB = LB;
cout << "-------------" << endl;    
cout << "--NCS ended--" << endl;
cout << "-------------" << endl;
cout << "Results: " << endl;    
cout << setprecision(12) << "UB: " << UB << endl;
cout << setprecision(12) << "LB: " << UB << endl;
cout << "time: " << setprecision(12) << tempo_decorrido(inicio_CPU) << endl;
return(0);
}

////////////////// run partial linear relaxation  ///////////////////
/////////////////////////////////////////////////////////////////////

cout << "Solve partial linear relaxation: " << endl;
for (int j = 0; j < m; j++)
{
  mod_denso.add(IloConversion(env_denso, y_denso[j], ILOBOOL));
}
cplex_denso.setParam(IloCplex::TiLim, tempo_limite_problema_estrategia_integralidade_parcial);
cplex_denso.setParam(IloCplex::Threads,1);
cplex_denso.setWarning(std::cerr);
cplex_denso.setOut(std::cout);
cplex_denso.solve();
cout << "time CP+PIS: " << tempo_decorrido(inicio_CPU) << endl;  
   for (int i = 0; i < n; i++)
{
  for (int j = 0; j < m; j++)
  {
      x_copia[i][j] = cplex_denso.getValue(x_denso[i][j])  ;
       
  }
  
}
// 
  for (int j = 0; j < m; j++)
  {
      y_copia[j] = cplex_denso.getValue(y_denso[j])  ;
  }
if (cplex_denso.getStatus() == IloAlgorithm::Optimal)
{
cout << "LB - PC + PIS: "<< setprecision(10) << cplex_denso.getObjValue() << endl;
LB = cplex_denso.getObjValue();

//se a solucao é inteira termina o algoritmo



solucao_inteira = true;
int conta_facilidades_parcialmente_abertas = 0;
for (int j = 0; j < m; j++)
{
    //cout << y_copia[j] << "(" << p[j] << ") " ;
    if ( (y_copia[j] > 0.00001) && (y_copia[j] < 0.99999) )
    {
           
           conta_facilidades_parcialmente_abertas = conta_facilidades_parcialmente_abertas + 1;
           //break;
    }
}
if (conta_facilidades_parcialmente_abertas > 0) 
{
    solucao_inteira = false;
}
else
{
    
   for (int i = 0; i < n; i++)
{
    for (int j = 0; j < m; j++)
    {
        if (  abs(round(x_copia[i][j]) - x_copia[i][j]) >= 0.0000001  )
        {
            solucao_inteira = false;
            break;
        }
    }
    if (solucao_inteira == false) break;
    
} 
}
if (solucao_inteira == true)
{
//algoritmo termina apos estrategia de integralidade parcial

// temrina algoritmo: solucao apos integralidade parcial é inteira

cout << "LB - CP+PIS : " << setprecision(10) <<  LB << endl; 
//fim_CPU = clock();
//double tempo_planos_cortantes = (double) (fim_CPU - inicio_CPU_pc_eip)/CLOCKS_PER_SEC;   
//fim algoritmo    
cout << "Integer Solution " << endl;
cout << endl;
solucao_inteira = true;

double otimo = LB;

}
}
else
{
 cout << "LB - PC + PIS: "<< setprecision(10) << cplex_denso.getBestObjValue() << endl;
LB = cplex_denso.getBestObjValue();   
}

fim_CPU = clock();

primeiro_lb = LB;
tempo_primeiro_lb = tempo_decorrido(inicio_CPU);

//exit(0);
if (solucao_inteira == true)
{
UB = LB;
cout << "-------------" << endl;    
cout << "--NCS ended--" << endl;
cout << "-------------" << endl;
cout << "Results: " << endl;    
cout << setprecision(12) << "UB: " << UB << endl;
cout << setprecision(12) << "LB: " << UB << endl;
cout << "time: " << setprecision(12) << tempo_decorrido(inicio_CPU) << endl;
return(0);

}

if (coeficientes_inteiros == true)
{
    
    if (  abs(ceil(LB) - LB) < 0.99999  )
    {
    LB = ceil(LB);
    }
      //cout << "LB apos arredondamento: "  << LB << endl;
     // primeiro_lb = LB;
    //tempo_primeiro_lb = tempo_decorrido(inicio_CPU);
}


//////////////////////////////////////////
//////////////////////////////////////////
//////////////////////////////////////////
////////// CUT - AND -  SOLVE PHASE  /////
//////////////////////////////////////////
//////////////////////////////////////////
//////////////////////////////////////////

double fo1 = 0;
double fo2 = 0;
for (int j = 0; j < m; j++)
{
   if (y_copia[j] > 0.99)
   {

   locais_promissores.push_back(j);
   soma_capacidades_locais_promissores = soma_capacidades_locais_promissores + p[j];

   }
   else
   {
       locais_nao_promissores.push_back(j); 
       //corte2[0].push_back(j);
       //fechadas.push_back(j);
   }
}
vector<int> vetor_conjunto_clientes;
vector<int> vetor_locais_promissores = cria_vetor_copia_de_lista(locais_promissores);

for (int i = 0; i < n; i++)
{
   vetor_conjunto_clientes.push_back(i);
}
cout << "_____________________________CUT-AND-SOLVE PHASE_________________________________________________" << endl;
bool solucao_inicial_heuristica_factivel = false;
while   (tempo_total < tempo_limite_total) 
{
quantidade_iteracoes = quantidade_iteracoes + 1;
if ( (coeficientes_inteiros == true) && (UB - ceil(LB) <= 0)  ) break;
if ( (coeficientes_inteiros == false) && ( UB - LB < 0.0001*LB )   ) break;
cout << "__________________________________________________________________________________________________________" << endl;    

cout << "Iter :" << quantidade_iteracoes << endl;
cout << endl;

fim_CPU = clock();
if (  tempo_decorrido(inicio_CPU) > tempo_limite_total )
{
cout << "-------------" << endl;    
cout << "--NCS ended--" << endl;
cout << "-------------" << endl;
cout << "Results: " << endl;    
cout << setprecision(12) << "UB: " << UB << endl;
cout << setprecision(12) << "LB: " << UB << endl;
cout << "time: " << setprecision(12) << tempo_decorrido(inicio_CPU) << endl;
return(0);
}

bool problema_esparso_infactivel = false;
roda_heuristica = true;

bool problema_esparso_e_gap = false;
verifica_se_problema_esparso_e_gap(locais_promissores, problema_esparso_e_gap, vetor_conjunto_clientes,  p, soma_demandas, vetor_locais_promissores, soma_capacidades_locais_promissores, pen);
//cout << "valor de iter: " << iter << endl;

if (quantidade_iteracoes > 1)
{
tempo_limite_problema_esparso = tempo_limite_total - tempo_decorrido(inicio_CPU);
SolveSparseProblem(problema_esparso_e_gap, solucao_inicial_heuristica_factivel, quantidade_iteracoes,  vetor_locais_promissores, f, vetor_conjunto_clientes,  c, p, pen, d, UB,  coeficientes_inteiros, tempo_limite_problema_esparso, iteracao_que_encontrou_otimo, solucao_heuristica, atende, coeficientes_cortes_x,  clientes_cortes_x, limitantes_cortes_x, variavel_y_associada_ao_corte,  coeficientes_cortes_y, facilidades_cortes_y, limitantes_cortes_y, melhor_atende, n, m, Qmin, lb_esparso, LB_modelo_otimo);
}
else
{
        if  ( locais_promissores.size() > U  )
        {
                                cout << "Run Heuristic: " << endl;
                                double soma_capacidades_solucao_inicial_heuristica = 0;
                                vector<list<int>> s_aux(m);

                                    for (int i = 0; i < n; i++)
                                    {
                                        atende[i] = 0;
                                        s_aux[0].push_back(i);
                                        for (int j = 1; j < m ; j++)
                                        {
                                        if   (  (x_copia[i][j] > 0) &&  (x_copia[i][j] >= x_copia[i][atende[i]]  ) )
                                        {
                                            if (  (x_copia[i][j] == x_copia[i][atende[i]] ) && (c[i][j] < c[i][atende[i]] )   )
                                                {
                                                    s_aux[atende[i]].remove(i);
                                                    atende[i] = j;
                                                    s_aux[j].push_back(i);
                                                }
                                            if  (x_copia[i][j] > x_copia[i][atende[i]] )
                                                {
                                                    s_aux[atende[i]].remove(i);
                                                    atende[i] = j;
                                                    s_aux[j].push_back(i);
                                                }   

                                        }
                                            
                                        }
                                    }
                                    
                                    for (int j = 0; j < m; j++)
                                    {
                                        if( (s_aux[j].empty()) && (esta_na_lista(locais_promissores, j)  )  )
                                        {
                                          //cout << "facilidade vazia e pertence aos locais promissores" << j << endl;
                                          // seleciona o cliente mais proximo de j:
                                          int mais_proximo = 0;
                                          for (int i = 1; i < n ; i++)
                                          {
                                            if (c[i][j] < c[mais_proximo][j]) mais_proximo = i;
                                          }
                                          atende[mais_proximo] = j;
                                                     
                                        }
                                    }
                                  
                                            abertas.clear();
                                            fechadas.clear();
                                            fo1 = 0;
                                            fo2 = 0;
                                                for (int j = 0; j < m; j++)
                                                {

                                                pen[j] = -p[j];
                                                s[j].clear();
                                                }
                                                    //atualiza fo
                                            for (int i = 0; i < n; i++)    
                                            {

                                            s[atende[i]].push_back(i);
                                            fo1 = fo1 + c[i][atende[i]]*d[i];
                                            pen[atende[i]] = pen[atende[i]] + d[i];

                                            }
                                                                            
                                            for (int j = 0; j < m; j++)
                                            {
                                                if (s[j].empty())
                                                {
                                                    fechadas.push_back(j);
                                                }
                                                else
                                                {
                                                    soma_capacidades_solucao_inicial_heuristica = soma_capacidades_solucao_inicial_heuristica + p[j];
                                                    abertas.push_back(j);
                                                    fo2 = fo2 + f[j];
                                                    if (pen[j] > 0)
                                                    {
                                                    factivel = false;
                                                    fo2 = fo2 + 100000*pen[j];
                                                //break;
                                                    }
                                                }
                                                
                                            }
                                        fo = fo1 + fo2; 
                                        cout << "Initial solution obj value: " << fo << endl;
                                        //cout << "facilidades abertas na solução da heurística:" << abertas.size() << endl;
                                        //cout << "soma capacidades das facilidades abertas na solução da heurística: " << soma_capacidades_solucao_inicial_heuristica << endl;
                                        //cout << "lista facilidades abertas apos a aplicacao da heuristica" << endl;
                                        //imprime_lista(abertas);
                                        //cout << endl;

                                bool solucao_inicial_heuristica_factivel = true;
                                bool resolveu_sp_na_otimalidade = false;
                                LocalSearch(m, n, p,  locais_promissores, pen,s,c,f,d,atende,coeficientes_inteiros,otimo,UB,abertas, fechadas, fo,soma_demandas, imprime_detalhes, solucao_inicial_heuristica_factivel , problema_esparso_e_gap, Qmin, R,  resolveu_sp_na_otimalidade, tempo_limite_heuristica, LB_modelo_otimo);
                                tempo_limite_problema_esparso = tempo_limite_total - tempo_decorrido(inicio_CPU);
                                cout << setprecision(12) << "Upper bound : " << fo << endl;
                                primeiro_UB = fo;
                                // tempo_primeiro_UB = tempo_decorrido(inicio_CPU);
                                // cout << "tempo decorrido: " << tempo_primeiro_UB << endl;
                                if (solucao_inicial_heuristica_factivel)
                                {
                                    primeiro_UB = fo;
                                   
                                    if(  (coeficientes_inteiros == false) && (primeiro_UB - LB < 0.0001*LB) )
                                    {
                                     UB = primeiro_UB;  
                                    cout << setprecision(12) << "UB: " << UB << endl;
                                    cout << setprecision(12) << "LB: " << UB << endl;
                                    cout << "time: " << setprecision(12) << tempo_decorrido(inicio_CPU) << endl;
                                    return(0);
                                                            


                            
                            return(0);

                                }
                            
                                    if( (coeficientes_inteiros == true) && (primeiro_UB - ceil(LB) <= 0) )
                                    {
                                   UB = primeiro_UB;
                                                                      cout << setprecision(12) << "UB: " << UB << endl;
                                    cout << setprecision(12) << "LB: " << UB << endl;
                                    cout << "time: " << setprecision(12) << tempo_decorrido(inicio_CPU) << endl;
                                    return(0);
                                }
                                }

                                                    if (fo < UB) 
                                                    {
                                                        UB = fo;
                                                        melhor_atende = atende;
                                                        iteracao_que_encontrou_otimo = quantidade_iteracoes; 
                                                    }
                    if (resolveu_sp_na_otimalidade == false) 
                    {
                        SolveSparseProblem(problema_esparso_e_gap, solucao_inicial_heuristica_factivel, quantidade_iteracoes,  vetor_locais_promissores, f, vetor_conjunto_clientes,  c, p, pen, d, UB,  coeficientes_inteiros, tempo_limite_problema_esparso, iteracao_que_encontrou_otimo, solucao_heuristica, atende, coeficientes_cortes_x,  clientes_cortes_x, limitantes_cortes_x, variavel_y_associada_ao_corte,  coeficientes_cortes_y, facilidades_cortes_y, limitantes_cortes_y, melhor_atende, n, m, U, lb_esparso, LB_modelo_otimo);

                    }
        }
        else
        {
          //  arq_utilizou_heuristica << 0 << endl; exit(0);
            tempo_limite_problema_esparso = tempo_limite_total - tempo_decorrido(inicio_CPU);
             SolveSparseProblem(problema_esparso_e_gap, solucao_inicial_heuristica_factivel, quantidade_iteracoes,  vetor_locais_promissores, f, vetor_conjunto_clientes,  c, p, pen, d, UB,  coeficientes_inteiros, tempo_limite_problema_esparso, iteracao_que_encontrou_otimo , solucao_heuristica, atende, coeficientes_cortes_x,  clientes_cortes_x, limitantes_cortes_x, variavel_y_associada_ao_corte,  coeficientes_cortes_y, facilidades_cortes_y, limitantes_cortes_y, melhor_atende, n, m, U, lb_esparso, LB_modelo_otimo);

                primeiro_UB = UB;
                tempo_primeiro_UB = tempo_decorrido(inicio_CPU);
                iteracao_que_encontrou_otimo = quantidade_iteracoes;
                 //arq_utilizou_heuristica << 0 << endl;
    
    
                    }


        }


if (coeficientes_inteiros == false)
{
if( UB - LB < 0.0001*LB)
{
cout << "-------------" << endl;    
cout << "--NCS ended--" << endl;
cout << "-------------" << endl;
cout << "Results: " << endl;    
cout << setprecision(12) << "UB: " << UB << endl;
cout << setprecision(12) << "LB: " << UB << endl;
cout << "time: " << setprecision(12) << tempo_decorrido(inicio_CPU) << endl;
return(0);
}
}
else
{
if( UB - ceil(LB) <= 0)
{
cout << "-------------" << endl;    
cout << "--NCS ended--" << endl;
cout << "-------------" << endl;
cout << "Results: " << endl;    
cout << setprecision(12) << "UB: " << UB << endl;
cout << setprecision(12) << "LB: " << UB << endl;
cout << "time: " << setprecision(12) << tempo_decorrido(inicio_CPU) << endl;
return(0);
}    
}
fim_CPU = clock();
if (locais_promissores.size() == m)
{
                                    cout << setprecision(12) << "UB: " << UB << endl;
                                    cout << setprecision(12) << "LB: " << UB << endl;
                                    cout << "time: " << setprecision(12) << tempo_decorrido(inicio_CPU) << endl;
                                    return(0);
    return(0);

}

fim_CPU = clock();
cout << setprecision(12) << "UB: " << UB << " time: " << (double) (fim_CPU - inicio_CPU)/CLOCKS_PER_SEC << endl;

cout << endl;
cout << "________________________________________________________________" << endl;
cout << "Solve RDP:" << endl;
cout << endl;
//inicio_CPU_auxiliar = clock();
resolve_denso:
    IloExpr raux (env_denso);
    for (std::list<int>::iterator k = locais_nao_promissores.begin() ; k !=  locais_nao_promissores.end() ; k++)
    {
    //cout << *k << " ";
     raux += y_denso[*k];
    }
    mod_denso.add( raux >= 1);
    raux.end();
cplex_denso.setParam(IloCplex::Threads,1);
cplex_denso.solve();

cout << "New LB: :" << cplex_denso.getObjValue() << endl; 
LB = cplex_denso.getObjValue();
if (coeficientes_inteiros == false)
{
if( UB - LB < 0.0001*LB)
{
cout << "-------------" << endl;    
cout << "--NCS ended--" << endl;
cout << "-------------" << endl;
cout << "Results: " << endl;    
cout << setprecision(12) << "UB: " << UB << endl;
cout << setprecision(12) << "LB: " << UB << endl;
cout << "time: " << setprecision(12) << tempo_decorrido(inicio_CPU) << endl;
return(0);
}
}
else
{
if( UB - ceil(LB) <= 0)
{
cout << "-------------" << endl;    
cout << "--NCS ended--" << endl;
cout << "-------------" << endl;
cout << "Results: " << endl;    
cout << setprecision(12) << "UB: " << UB << endl;
cout << setprecision(12) << "LB: " << UB << endl;
cout << "time: " << setprecision(12) << tempo_decorrido(inicio_CPU) << endl;
return(0);
}    
}
vector<int> facilidades_cheias;
vector<int> facilidades_fracionarias;
vector<int> facilidades_fechadas;
double soma_capacidades_cheias = 0; 
for (int j = 0; j < m; j++)
{
    y_copia[j] = cplex_denso.getValue(y_denso[j]);
      
      if (y_copia[j] > 0.00001  )
      { 
          if (y_copia[j] < 0.99999)
          {
            facilidades_fracionarias.push_back(j);
          }
          else
          {
             facilidades_cheias.push_back(j);
             soma_capacidades_cheias = soma_capacidades_cheias + p[j];
          }
            
      }  
      else
      {
         facilidades_fechadas.push_back(j);
      }
       
    for (int i = 0; i < n; i++)
    {
      x_copia[i][j] = cplex_denso.getValue(x_denso[i][j])  ;
    }
  
}
iter = iter + 1;
locais_promissores.clear();
locais_nao_promissores.clear();
vetor_locais_promissores.clear();

for (int j = 0; j < m; j++)
{
   if (cplex_denso.getValue(y_denso[j]) > 0.99)
   {
   locais_promissores.push_back(j);
   }
   else
   {
       locais_nao_promissores.push_back(j);
   }
}
vetor_locais_promissores = cria_vetor_copia_de_lista(locais_promissores);
bool solucao_denso_inteira = true;
for (int i = 0; i < n; i++)
{
    for (int j = 0; j < m; j++)
    {
       if (   abs(   cplex_denso.getValue(x_denso[i][j]) - round(cplex_denso.getValue(x_denso[i][j]))  ) > 0.000001 )
        {
            solucao_denso_inteira = false;
            goto label;
        }
    }
    
}
label:
if ( (solucao_denso_inteira == true)  && (LB <= UB))
{
        UB = LB;
cout << "-------------" << endl;    
cout << "--NCS ended--" << endl;
cout << "-------------" << endl;
cout << "Results: " << endl;    
cout << setprecision(12) << "UB: " << UB << endl;
cout << setprecision(12) << "LB: " << UB << endl;
cout << "time: " << setprecision(12) << tempo_decorrido(inicio_CPU) << endl;
return(0);
}
if (primeira_solucao_factivel_encontrada == false)
{
// atualiza x_copia e y_copia
for (int j = 0; j < m; j++)
{
    y_copia[j] = cplex_denso.getValue(y_denso[j]);
    for (int i = 0; i < n; i++)
    {
          x_copia[i][j] = cplex_denso.getValue(x_denso[i][j]);
    }
    
}
}


LB = cplex_denso.getObjValue();
if (coeficientes_inteiros == true)
{
    LB = ceil(LB);
}

cout << "-----------------Iteration finished-------------------" << endl;
cout << setprecision(12) << "New LB: " << LB << endl; //alterado
cout << setprecision(12) << "current UB: " << UB << endl;
cout << setprecision(12) << "time: " << tempo_decorrido(inicio_CPU) << endl;
cout << "------------------------------------------------------" << endl;

}

if (coeficientes_inteiros == false)
{
if( UB - LB < 0.0001*LB)
{
cout << "-------------" << endl;    
cout << "--NCS ended--" << endl;
cout << "-------------" << endl;
cout << "Results: " << endl;    
cout << setprecision(12) << "UB: " << UB << endl;
cout << setprecision(12) << "LB: " << UB << endl;
cout << "time: " << setprecision(12) << tempo_decorrido(inicio_CPU) << endl;
return(0);
}
}
else
{
if( UB - ceil(LB) <= 0)
{
cout << "-------------" << endl;    
cout << "--NCS ended--" << endl;
cout << "-------------" << endl;
cout << "Results: " << endl;    
cout << setprecision(12) << "UB: " << UB << endl;
cout << setprecision(12) << "LB: " << UB << endl;
cout << "time: " << setprecision(12) << tempo_decorrido(inicio_CPU) << endl;
return(0);
}    
}
if (tempo_total > tempo_limite_total)
{
cout << "-------------" << endl;    
cout << "--NCS ended--" << endl;
cout << "-------------" << endl;
cout << "Results: " << endl;    
cout << setprecision(12) << "UB: " << UB << endl;
cout << setprecision(12) << "LB: " << UB << endl;
cout << "time: " << setprecision(12) << tempo_decorrido(inicio_CPU) << endl;
return(0);
}

		}
	catch (IloException& ex) {
		cerr << "Error: " << ex << endl;
	}
}//fim algoritmo
