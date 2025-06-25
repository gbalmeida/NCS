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
#include <ilcplex/ilocplex.h>

using namespace std;
ILOSTLBEGIN



void CuttingPlane(int m, int n, vector<double> f, vector<double> p, vector<vector<double>> c, double &LB, double soma_demandas,  vector<double> d,  vector<double> &y_copia, vector<vector<double>>  &x_copia, double capacidade_total, double demanda_total, clock_t &fim_CPU, clock_t &inicio_CPU , int &quantidade_iteracoes, int &iteracao_que_encontrou_otimo,  double &tempo_resolve_problema_denso,  double tolerancia_pc,  int bigM, bool aplica_funcoes_de_lifting, bool aplica_lpi, bool aplica_lci_geral, bool aplica_fenchel_cuts,  list<vector<double>> &coeficientes_cortes_x, list<vector<int>> &clientes_cortes_x, vector<double> &limitantes_cortes_x, vector<int> &variavel_y_associada_ao_corte, list<vector<double>> &coeficientes_cortes_y, list<vector<int>> &facilidades_cortes_y, vector<double> &limitantes_cortes_y,  int max_iter_fc,  bool aplica_heuristica_gulosa1, bool aplica_heuristica_gulosa2, bool aplica_separacao_exata_de_cis, bool aplica_vubs, bool CIgreedyNew, int tempo_maximo_pc,  double tempo_limite_problema_estrategia_integralidade_parcial,   IloEnv& env_denso, IloModel& mod_denso, IloCplex& cplex_denso, IloNumVarArray& y_denso, IloArray<IloNumVarArray> x_denso, bool &solucao_CPPIS_inteira, IloConstraintArray &restricoes_fcs)
{

bool imprime_logs_cplex = false; 
clock_t inicio = clock();
double tempo_resolve_relaxacao = 0;
clock_t inicio_CPU_pc_eip = clock();
IloExpr expfo_denso(env_denso);
for (int j = 0; j <  m; j++){
	expfo_denso += f[j] * y_denso[j];
	for(int i = 0; i < n; i++){
		expfo_denso += c[i][j] * x_denso[i][j] * d[i];
	}
}  
IloAdd(mod_denso, IloMinimize(env_denso, expfo_denso));
expfo_denso.end();
// assigment constraints  *****************************************************************
for (int i = 0; i < n; i++){
			IloExpr r22(env_denso);
			for (int j = 0; j < m; j++){
				r22 += x_denso[i][j];
                if ( d[i] > p[j] )
                {
                 mod_denso.add( x_denso[i][j] <= 0)   ;
                }

			}
			mod_denso.add(r22 == 1);
			r22.end();
		} 
// capacity constraints *****************************************************************
			   for (int j = 0; j < m; j++){
			IloExpr r32(env_denso);
			for (int i = 0; i < n; i++){
				r32 += d[i]*x_denso[i][j];
			}
			mod_denso.add(r32 <= p[j]*y_denso[j]);
			r32.end();
		} 
// demanda total *****************************************************************
IloExpr r34(env_denso);
			   for (int j = 0; j < m; j++){
			
				r34 += p[j]*y_denso[j];
			}
			mod_denso.add(r34 >= demanda_total);
			r34.end();
            IloExpr r35(env_denso);

cplex_denso.setWarning(env_denso.getNullStream()); // Eliminar warnings
cplex_denso.setOut(env_denso.getNullStream()); /// Eliminar os logs do solver
cplex_denso.setParam(IloCplex::Threads,1);
cplex_denso.setParam(IloCplex::Param::Parallel,1);

cout << "Running Cutting-Plane method " << endl;
	
cplex_denso.solve();
clock_t fim = clock();
tempo_resolve_relaxacao = tempo_resolve_relaxacao + (double) (fim - inicio)/CLOCKS_PER_SEC;
cout << "LB - Linear Relaxation:" << setprecision(12) << cplex_denso.getObjValue() << endl;  

IloConstraintArray restricoes_vubs(env_denso);
IloConstraintArray restricoes_lcis(env_denso);
//IloConstraintArray restricoes_fcs(env_denso);

vector<int> Itens;
for (int i = 0; i < n; i++)
{
    Itens.push_back(i);
}

int total_tentativas = 0;
int tentativas_opcao0 = 0;

vector<int> facilidades_cheias;
vector<int> facilidades_fracionarias;
vector<int> facilidades_fechadas;
double soma_capacidades_cheias = 0; 
double media_valores_fracionarias = 0;
for (int j = 0; j < m; j++)
{
   
    y_copia[j] = cplex_denso.getValue(y_denso[j]);
      
      if (y_copia[j] > 0.000001  )
      { 
          if (y_copia[j] < 0.999999)
          {
            facilidades_fracionarias.push_back(j);
          }
          else
          {
             facilidades_cheias.push_back(j);
             soma_capacidades_cheias = soma_capacidades_cheias + p[j];
          }
          media_valores_fracionarias = media_valores_fracionarias + y_copia[j];
            
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

    vector<double> y_copia_antes(m);
    vector<vector<double>> x_copia_antes(n, vector<double>(m) )   ;
list<int> locais_tipo1; 
bool pare_planos_cortantes_mais_eip = false;
int conta_iter = 0;
bool converteu_variaveis = false;
double tempo_planos_cortantes_mais_parcial = 0;
double melhor_LB = cplex_denso.getObjValue();
double LB_antes = cplex_denso.getObjValue();
bool y_copia_e_igual = false;
vector<int> solucao(n);

LB_antes = melhor_LB;
bool pare_planos_cortantes = false;
double li_antes;
double li_depois;
while (pare_planos_cortantes == false)
{

vector<double> y_copia_antes = y_copia;
vector<vector<double>> x_copia_antes = x_copia ;
li_antes = cplex_denso.getObjValue();
li_depois = li_antes;
bool roda_cplex = false;
int conta_cortes_efetivos = 0;
locais_tipo1.clear();
clock_t tempo_inicio = clock();

if (aplica_vubs == true)
{
for (int j = 0; j < m; j++)
   {

for (int i = 0; i < n; i++)
{
      if ( x_copia[i][j]  > y_copia[j] )
      {

      IloConstraint cons = y_denso[j] >= x_denso[i][j];
      mod_denso.add(cons); 
      restricoes_vubs.add(cons);
      conta_cortes_efetivos = conta_cortes_efetivos + 1;
      }
   
   }
}
}
clock_t tempo_fim = clock();
//////////////////////////////////
//////////////////////////////////
////////////////////////////////// APLICA CORTES SOBRE AS VARIÁVEIS X 
//////////////////////////////////
///////////////////////////////////


//seleciona locais tipo 1, atualiza y copia e x copia e aplica cortes de fluxo simples 
tempo_inicio = clock();
for (int j = 0; j < m; j++)
{
              if ( (y_copia[j] > 0.00001) && (y_copia[j] < 0.99999)  )
      {    
          locais_tipo1.push_back(j);
      }
      else
      {
            if (y_copia[j] > 0.00001)
            {
            for (int i = 0; i < n; i++)
                     {
              if ( (x_copia[i][j] > 0.000001) && (x_copia[i][j] < 0.999999)  ) //fonte unica violada. 
                {
                locais_tipo1.push_back(j);
                break;
                }
                }
            }

      }
}

if (locais_tipo1.size() > 0)
{
locais_tipo1.sort();
locais_tipo1.unique();
vector<int> vetor_locais_tipo1 = cria_vetor_copia_de_lista(locais_tipo1);
for (int j = 0; j < vetor_locais_tipo1.size(); j++)
{
double y_til = y_copia[vetor_locais_tipo1[j]];
vector<int> clientes_atendidos; 
bool frac = false;
if ( abs(round(y_til) - y_til) > 0.00001 )
{
    frac = true;
}
double soma_demandas_clientes_atendidos = 0;
int demanda_que_chega = 0;
double soma_demandas_da_cobertura = 0;
double soma_demandas_da_cobertura_minimal = 0;
double soma_demandas_pack = 0;
double capacidade_residual = p[vetor_locais_tipo1[j]];
list<int> cobertura;
list<int> cobertura_minimal;
list<int> pack;
double soma = 0;
bool violada = false;
vector<int> vetor_cobertura;
vector<int> vetor_cobertura_minimal;
vector<int> vetor_pack;
vector<int> vetor_D;
list<int> lista_D;
vector<int> vetor_N1;
vector<int> vetor_N2;
double soma_demandas_D = 0;
vector<double> x(n);
    for (int i = 0; i < n; i++)
    {
        x[i] = x_copia[i][vetor_locais_tipo1[j]];
    }
        double B = p[vetor_locais_tipo1[j]];
    vector<int> clientes_corte;
    vector<double> coef;
    double limitante = 0;
    double limitante_minimal = 0;
vector<double> aux_denso(n);
vector<double> aux3(n);
for (int i = 0; i < n; i++)
{
    aux_denso[i] = x_copia[i][vetor_locais_tipo1[j]];
    aux3[i] = (double) (1-x_copia[i][vetor_locais_tipo1[j]])/d[i];
    if (x_copia[i][vetor_locais_tipo1[j]] > 0.000001)
    {
        capacidade_residual = capacidade_residual - x_copia[i][vetor_locais_tipo1[j]]*d[i] ;
        clientes_atendidos.push_back(i);
        soma_demandas_clientes_atendidos = soma_demandas_clientes_atendidos + d[i];
        demanda_que_chega = demanda_que_chega + d[i]*x_copia[i][vetor_locais_tipo1[j]];
        if (x_copia[i][vetor_locais_tipo1[j]] > 0.999999)
        {
            vetor_D.push_back(i);
            lista_D.push_back(i);
            soma_demandas_D = soma_demandas_D + d[i];
        }
        else
        {
             vetor_N1.push_back(i);
        }
     }
    else
    {
        vetor_N2.push_back(i);
    }
    
}
funcao_ordena_decrescente2(vetor_N1, x, d);
funcao_ordena_decrescente(vetor_N2, d);
vector<int> best_CI;
double  valor_best_CI = 10e10;
double soma_demandas_best_CI = 0;
vetor_cobertura = vetor_D;
soma_demandas_da_cobertura = soma_demandas_D;  
vector<int> vetor_cobertura_auxiliar = vetor_D;
vector<int> vetor_pack_auxiliar = vetor_D;
list<int> cobertura_auxiliar = lista_D;
list<int> pack_auxiliar = lista_D;
double soma_demandas_da_cobertura_auxiliar = soma_demandas_D;
double soma_demandas_pack_auxiliar = soma_demandas_D;
bool corte_de_cobertura_encontrado = false;
bool cobertura_encontrada = false;
double capacidade_residual_D = p[vetor_locais_tipo1[j]] - soma_demandas_D;
bool tentativa_1 = true;
bool tentativa_2 = false;
bool tentativa_1_deu_certo = false;
bool tentativa_1_pack_deu_certo = false;
bool movimento_acrescenta1_cliente = false;
bool movimento_acrescenta2_cliente = false;
bool movimento_acrescenta3_cliente = false;
tempo_inicio = clock();
if ( (aplica_heuristica_gulosa1 == true) && (violada == false))
{
CIgreedyA( vetor_N1, vetor_N2, aux_denso, soma_demandas_da_cobertura, B, d, vetor_cobertura, corte_de_cobertura_encontrado, x, y_til, best_CI, valor_best_CI);

double soma_aux = 0;
for (int i = 0; i < vetor_cobertura.size(); i++)
{
    cobertura.push_back(vetor_cobertura[i]);
    
    if (soma_aux + d[vetor_cobertura[i]] < B )
    {
        pack.push_back(vetor_cobertura[i]);
        vetor_pack.push_back(vetor_cobertura[i]);
        soma_demandas_pack = soma_demandas_pack + d[vetor_cobertura[i]];
        soma_aux = soma_aux + d[vetor_cobertura[i]];
    }
    
}
limitante = (vetor_cobertura.size() - 1);
 cobertura_encontrada = true;
if (corte_de_cobertura_encontrado == false)
{
    tentativa_2 = true;
}
}

if (aplica_heuristica_gulosa1 == false) tentativa_2 == true; 
if ( (tentativa_2 == true) )
{
vetor_cobertura_auxiliar = vetor_D;
soma_demandas_da_cobertura_auxiliar = soma_demandas_D;

if (aplica_heuristica_gulosa2 == true)
{
CIgreedyB( vetor_N1, vetor_N2, aux_denso, soma_demandas_da_cobertura_auxiliar, B, d, vetor_cobertura_auxiliar, corte_de_cobertura_encontrado, x, y_til, best_CI, valor_best_CI);
if (corte_de_cobertura_encontrado == true)
{
  
    cobertura_encontrada = true;
      vetor_cobertura = vetor_cobertura_auxiliar;
    soma_demandas_da_cobertura = soma_demandas_da_cobertura_auxiliar;
    double soma_aux = 0;
    pack.clear();
    cobertura.clear();
    vetor_pack.clear();
    soma_demandas_pack = 0;
    for (int i = 0; i < vetor_cobertura.size(); i++)
    {
    cobertura.push_back(vetor_cobertura[i]);
        if (soma_aux + d[vetor_cobertura[i]] < B )
    {
        pack.push_back(vetor_cobertura[i]);
        vetor_pack.push_back(vetor_cobertura[i]);
        soma_demandas_pack = soma_demandas_pack + d[vetor_cobertura[i]];
        soma_aux = soma_aux + d[vetor_cobertura[i]];
    }
        }
    limitante = vetor_cobertura.size() - 1;
    }
else
{
 movimento_acrescenta1_cliente = true;
 if (aplica_heuristica_gulosa1 == false)
 {
       vetor_cobertura = vetor_cobertura_auxiliar;
    soma_demandas_da_cobertura = soma_demandas_da_cobertura_auxiliar;
    double soma_aux = 0;
    pack.clear();
    cobertura.clear();
    vetor_pack.clear();
    soma_demandas_pack = 0;
    for (int i = 0; i < vetor_cobertura.size(); i++)
    {
    cobertura.push_back(vetor_cobertura[i]);
   
    if (soma_aux + d[vetor_cobertura[i]] < B )
    {
        pack.push_back(vetor_cobertura[i]);
        vetor_pack.push_back(vetor_cobertura[i]);
        soma_demandas_pack = soma_demandas_pack + d[vetor_cobertura[i]];
        soma_aux = soma_aux + d[vetor_cobertura[i]];
    }
    
    }
    limitante = vetor_cobertura.size() - 1;
 }
}
}
else
{
    movimento_acrescenta1_cliente = true;
}
}
if (CIgreedyNew == false)
{
    movimento_acrescenta1_cliente = false;
    movimento_acrescenta2_cliente = false;
    movimento_acrescenta3_cliente = false;
}
tempo_inicio = clock();
if (  (movimento_acrescenta1_cliente == true) && (vetor_N1.size() >=1   ) )
{
vetor_cobertura_auxiliar = vetor_D;
soma_demandas_da_cobertura_auxiliar = soma_demandas_D;
if (CIgreedyNew == true)
{
AddOneCostumer( vetor_N1, vetor_N2, aux_denso, soma_demandas_da_cobertura_auxiliar, B, d, vetor_cobertura_auxiliar, corte_de_cobertura_encontrado, x, y_til, best_CI, valor_best_CI);
}
if (corte_de_cobertura_encontrado == true)
{
   // cout << "acrescenta 1 cliente deu certo" << endl;
    cobertura_encontrada = true;
      vetor_cobertura = vetor_cobertura_auxiliar;
    soma_demandas_da_cobertura = soma_demandas_da_cobertura_auxiliar;
    double soma_aux = 0;
    pack.clear();
    cobertura.clear();
    vetor_pack.clear();
    soma_demandas_pack = 0;
    for (int i = 0; i < vetor_cobertura.size(); i++)
    {
    cobertura.push_back(vetor_cobertura[i]);
    if (soma_aux + d[vetor_cobertura[i]] < B )
    {
        pack.push_back(vetor_cobertura[i]);
        vetor_pack.push_back(vetor_cobertura[i]);
        soma_demandas_pack = soma_demandas_pack + d[vetor_cobertura[i]];
        soma_aux = soma_aux + d[vetor_cobertura[i]];
    }
    }
    limitante = vetor_cobertura.size() - 1;
}
else
{
 movimento_acrescenta2_cliente = true;
}
}
if (  (movimento_acrescenta2_cliente == true) && (vetor_N1.size() >=2   ) )
{
vetor_cobertura_auxiliar = vetor_D;
soma_demandas_da_cobertura_auxiliar = soma_demandas_D;
if (CIgreedyNew == true)
{
AddTwoCostumers( vetor_N1, vetor_N2, aux_denso, soma_demandas_da_cobertura_auxiliar, B, d, vetor_cobertura_auxiliar, corte_de_cobertura_encontrado, x, y_til, best_CI, valor_best_CI);
}
if (corte_de_cobertura_encontrado == true)
{
     cobertura_encontrada = true;
    vetor_cobertura = vetor_cobertura_auxiliar;
    soma_demandas_da_cobertura = soma_demandas_da_cobertura_auxiliar;
    double soma_aux = 0;
    pack.clear();
    cobertura.clear();
    vetor_pack.clear();
    soma_demandas_pack = 0;
    for (int i = 0; i < vetor_cobertura.size(); i++)
    {
    cobertura.push_back(vetor_cobertura[i]);
    if (soma_aux + d[vetor_cobertura[i]] < B )
    {
        pack.push_back(vetor_cobertura[i]);
        vetor_pack.push_back(vetor_cobertura[i]);
        soma_demandas_pack = soma_demandas_pack + d[vetor_cobertura[i]];
        soma_aux = soma_aux + d[vetor_cobertura[i]];
    }
    }
    limitante = vetor_cobertura.size() - 1;
   }
else
{
 movimento_acrescenta3_cliente = true;
}
}
if (  (movimento_acrescenta3_cliente == true) && (vetor_N1.size() >=3   ) )
{
vetor_cobertura_auxiliar = vetor_D;
soma_demandas_da_cobertura_auxiliar = soma_demandas_D;
if (CIgreedyNew == true)
{
AddThreeCostumers( vetor_N1, vetor_N2, aux_denso, soma_demandas_da_cobertura_auxiliar, B, d, vetor_cobertura_auxiliar, corte_de_cobertura_encontrado, x, y_til, best_CI, valor_best_CI);
}
if (corte_de_cobertura_encontrado == true)
{
     cobertura_encontrada = true;
    vetor_cobertura = vetor_cobertura_auxiliar;
    soma_demandas_da_cobertura = soma_demandas_da_cobertura_auxiliar;
    double soma_aux = 0;
    pack.clear();
    cobertura.clear();
    vetor_pack.clear();
    soma_demandas_pack = 0;
    for (int i = 0; i < vetor_cobertura.size(); i++)
    {
    cobertura.push_back(vetor_cobertura[i]);
    
    if (soma_aux + d[vetor_cobertura[i]] < B )
    {
        pack.push_back(vetor_cobertura[i]);
        vetor_pack.push_back(vetor_cobertura[i]);
        soma_demandas_pack = soma_demandas_pack + d[vetor_cobertura[i]];
        soma_aux = soma_aux + d[vetor_cobertura[i]];
    }
    }
    limitante = vetor_cobertura.size() - 1;
}
}
if ( (aplica_separacao_exata_de_cis == true) && (corte_de_cobertura_encontrado == false)   )
{

vetor_cobertura_auxiliar = vetor_D;
soma_demandas_da_cobertura_auxiliar = soma_demandas_D;
cobertura_encontrada = false;
if (aplica_separacao_exata_de_cis == true)
{
CIexact( B, vetor_N1,  vetor_N2, vetor_D , y_til, x, d, vetor_cobertura_auxiliar, soma_demandas_da_cobertura_auxiliar, corte_de_cobertura_encontrado, soma_demandas_D, cobertura_encontrada, best_CI, valor_best_CI);
}

 if (corte_de_cobertura_encontrado == true)
 {  
    cobertura_encontrada = true;
     vetor_cobertura = vetor_cobertura_auxiliar;
    soma_demandas_da_cobertura = soma_demandas_da_cobertura_auxiliar;
    double soma_aux = 0;
    pack.clear();
    cobertura.clear();
    vetor_pack.clear();
    soma_demandas_pack = 0;
    for (int i = 0; i < vetor_cobertura.size(); i++)
    {
    cobertura.push_back(vetor_cobertura[i]);
    if (soma_aux + d[vetor_cobertura[i]] < B )
    {
        pack.push_back(vetor_cobertura[i]);
        vetor_pack.push_back(vetor_cobertura[i]);
        soma_demandas_pack = soma_demandas_pack + d[vetor_cobertura[i]];
        soma_aux = soma_aux + d[vetor_cobertura[i]];
    }
    }
    limitante = vetor_cobertura.size() - 1;
    double aux2 = 0;
    for (int i = 0; i < vetor_cobertura_auxiliar.size(); i++)
    {
        aux2 = aux2 + x[vetor_cobertura_auxiliar[i]];
    }
 }

else
{
    if (cobertura_encontrada == true)
    {
 
    vetor_cobertura = vetor_cobertura_auxiliar;
    soma_demandas_da_cobertura = soma_demandas_da_cobertura_auxiliar;
    double soma_aux = 0;
    pack.clear();
    cobertura.clear();
    vetor_pack.clear();
    soma_demandas_pack = 0;
    for (int i = 0; i < vetor_cobertura.size(); i++)
    {
    cobertura.push_back(vetor_cobertura[i]);
    if (soma_aux + d[vetor_cobertura[i]] < B )
    {
        pack.push_back(vetor_cobertura[i]);
        vetor_pack.push_back(vetor_cobertura[i]);
        soma_demandas_pack = soma_demandas_pack + d[vetor_cobertura[i]];
        soma_aux = soma_aux + d[vetor_cobertura[i]];
    }

    }
    limitante = vetor_cobertura.size() - 1;
  
    double aux2 = 0;
    for (int i = 0; i < vetor_cobertura_auxiliar.size(); i++)
    {
        aux2 = aux2 + x[vetor_cobertura_auxiliar[i]];
    }




    }
    else  //não vai usar CI ou LCI. Constroi o pack
    {

    double soma_aux = 0;
    bool pack_construido = false;
    //cobertura.clear();
    vetor_pack = vetor_D;
    soma_demandas_pack = soma_demandas_D;
    soma_aux = soma_demandas_D;
    for (int i = 0; i < vetor_N1.size(); i++)
    {
 
    if (soma_aux + d[vetor_N1[i]] < B )
    {
        //pack.push_back(vetor_D[i]);
        vetor_pack.push_back(vetor_N1[i]);
        soma_demandas_pack = soma_demandas_pack + d[vetor_N1[i]];
        soma_aux = soma_aux + d[vetor_N1[i]];
    }
    else
    {
        pack_construido = true;
        break;
    }

    }

    if (pack_construido == false)
    {
    for (int i = 0; i < vetor_N2.size(); i++)
    {
        
    
    if (soma_aux + d[vetor_N2[i]] < B )
    {
        //pack.push_back(vetor_D[i]);
        vetor_pack.push_back(vetor_N2[i]);
        soma_demandas_pack = soma_demandas_pack + d[vetor_N2[i]];
        soma_aux = soma_aux + d[vetor_N2[i]];
    }
    else
    {
        pack_construido = true;
        break;
    }

    }
    }
   
        for (int i = 0; i < vetor_pack.size(); i++)
    {
        pack.push_back(vetor_pack[i]);
    }

    

    }
    

}

}
if (   (cobertura_encontrada == true) && (aplica_funcoes_de_lifting == true) && (violada == false) /*&& (soma_demandas_clientes_atendidos > p[vetor_locais_tipo1[j]])*/  )
{

clock_t tempo_inicio = clock();
clientes_corte.clear();
coef.clear();
    bool corte_violado = false;
    LCIsuperaditivefunction(d , Itens, B, vetor_cobertura, coef, clientes_corte, limitante, x, corte_violado, y_til);
    double soma_aux = 0;
    for (int i = 0; i < coef.size(); i++)
    {
        soma_aux = soma_aux + x[clientes_corte[i]]*coef[i];
    }
if (corte_violado == true)
{
    violada = true;
    IloExpr rx(env_denso);
        clientes_cortes_x.push_back(clientes_corte);
        coeficientes_cortes_x.push_back(coef);
        limitantes_cortes_x.push_back(limitante);
        variavel_y_associada_ao_corte.push_back(vetor_locais_tipo1[j]); 
    for (int i = 0; i < clientes_corte.size(); i++)
    {
        rx += x_denso[clientes_corte[i]][vetor_locais_tipo1[j]]*coef[i];
    }
    IloConstraint cons = rx <= limitante*y_denso[vetor_locais_tipo1[j]];
    mod_denso.add(cons); 
    restricoes_lcis.add(cons);
    rx.end();
    conta_cortes_efetivos = conta_cortes_efetivos + 1;
}
clock_t tempo_fim = clock();
}// fim lifting lethford
if (   (cobertura_encontrada == true) &&   (aplica_lci_geral == true) && (violada == false) /*&& (soma_demandas_clientes_atendidos > p[vetor_locais_tipo1[j]])*/   )
{
clock_t tempo_inicio = clock();
clientes_corte.clear();
coef.clear();
   LCIsequentiallifting(d , Itens,B, vetor_cobertura, coef,clientes_corte, limitante, x, y_til);
if (verifica_corte_e_violado(d ,B, coef, clientes_corte,limitante, x, y_til) )
{
        violada = true;
        IloExpr rx(env_denso);
        clientes_cortes_x.push_back(clientes_corte);
        coeficientes_cortes_x.push_back(coef);
        limitantes_cortes_x.push_back(limitante);
        variavel_y_associada_ao_corte.push_back(vetor_locais_tipo1[j]);   
        for (int i = 0; i < clientes_corte.size(); i++)
        {
            rx += x_denso[clientes_corte[i]][vetor_locais_tipo1[j]]*coef[i];
        }
            IloConstraint cons = rx <= limitante*y_denso[vetor_locais_tipo1[j]];
            mod_denso.add(cons); 
            restricoes_lcis.add(cons);
        rx.end();
              conta_cortes_efetivos = conta_cortes_efetivos + 1;
}// fim lci geral
else
{
if (corte_de_cobertura_encontrado == true)
{
    violada = true;
    clientes_corte.clear();
    clientes_corte = vetor_cobertura;
    coef.clear();
    for (int i = 0; i < vetor_cobertura.size(); i++)
    {
        coef.push_back(1);
    }
    limitante = vetor_cobertura.size() - 1;
        clientes_cortes_x.push_back(clientes_corte);
        coeficientes_cortes_x.push_back(coef);
        limitantes_cortes_x.push_back(limitante);
        variavel_y_associada_ao_corte.push_back(vetor_locais_tipo1[j]); 
            IloExpr rx(env_denso);
    for (int i = 0; i < clientes_corte.size(); i++)
    {
        rx += x_denso[clientes_corte[i]][vetor_locais_tipo1[j]]*coef[i];
    }
            IloConstraint cons = rx <= limitante*y_denso[vetor_locais_tipo1[j]];
            mod_denso.add(cons); 
            restricoes_lcis.add(cons);
    rx.end();
          conta_cortes_efetivos = conta_cortes_efetivos + 1;
          }
 
}
}
if (   (aplica_lpi == true)  && (violada == false) /*&& (violada == false) && (soma_demandas_clientes_atendidos > p[vetor_locais_tipo1[j]]) */)
{
vector<int> vetor_clientes_restantes;
vector<double> alfa;
vector<double> lambda;
bool corte_violado = false;
LPIsuperaditivefunction(d , Itens, B, vetor_pack, pack, alfa, lambda, vetor_clientes_restantes, limitante, x, corte_violado, soma_demandas_pack);
soma = 0;
vector<int> clientes_corte;
vector<double> coef;
double soma_lambdas;
for (int i = 0; i < vetor_clientes_restantes.size(); i++)
{
    if (  abs(alfa[i]) > 0.000001)
    {
    clientes_corte.push_back(vetor_clientes_restantes[i]);
    coef.push_back(alfa[i]);
    }
}
for (int i = 0; i < vetor_pack.size(); i++)
{
    if (   abs(lambda[i]) > 0.0000001 )
    {
    clientes_corte.push_back(vetor_pack[i]);
    coef.push_back(-lambda[i]);
    soma_lambdas = soma_lambdas + lambda[i];
    }
}
double limitante_lpi = -soma_lambdas;
for (int i = 0; i < clientes_corte.size(); i++)
{
    soma = soma + x[clientes_corte[i]]*coef[i];
}
if (soma > limitante_lpi*y_til)
{
      violada = true;
        clientes_cortes_x.push_back(clientes_corte);
        coeficientes_cortes_x.push_back(coef);
        limitantes_cortes_x.push_back(limitante_lpi);
        variavel_y_associada_ao_corte.push_back(vetor_locais_tipo1[j]); 
    IloExpr rx(env_denso);
    for (int i = 0; i < clientes_corte.size(); i++)
    {
        rx += x_denso[clientes_corte[i]][vetor_locais_tipo1[j]]*coef[i];
       
    }
    mod_denso.add(rx <= limitante_lpi*y_denso[vetor_locais_tipo1[j]]);
    rx.end();
             conta_cortes_efetivos = conta_cortes_efetivos + 1;
}
}
if ( (violada == false) && (y_til > 0.999999) && (aplica_fenchel_cuts == true )  )
{
clock_t tempo_inicio = clock();
bool impossivel_cortar = true;
clientes_corte.clear();
coef.clear();
double limitante = 0;
bool corte_violado = false;
FCColumnGeneration(d , Itens,  B, coef, clientes_corte, limitante, x, corte_violado, impossivel_cortar, tentativas_opcao0, total_tentativas,  bigM, max_iter_fc);
if (corte_violado == true)
{
    violada = true;
    IloExpr rx(env_denso);
                //guarda o corte
        clientes_cortes_x.push_back(clientes_corte);
        coeficientes_cortes_x.push_back(coef);
        limitantes_cortes_x.push_back(limitante);
        variavel_y_associada_ao_corte.push_back(vetor_locais_tipo1[j]); 
    for (int i = 0; i < clientes_corte.size(); i++)
    {
        rx += x_denso[clientes_corte[i]][vetor_locais_tipo1[j]]*coef[i];
       
    }
                IloConstraint cons = rx <= limitante*y_denso[vetor_locais_tipo1[j]];
            mod_denso.add(cons); 
            restricoes_fcs.add(cons);
    rx.end();
          conta_cortes_efetivos = conta_cortes_efetivos + 1;
}
} //fim para cada vetor tipo 1
}
} //fim  "se há pelo menos um local tipo 1"
bool solucao_inteira = true;
if ( (conta_cortes_efetivos > 0   ) /*||  (roda_cplex == true)*/ )
{
cplex_denso.setParam(IloCplex::Threads,1);
tempo_inicio = clock();
cplex_denso.solve();
tempo_fim = clock();
tempo_resolve_relaxacao = tempo_resolve_relaxacao + (double) (tempo_fim - tempo_inicio)/CLOCKS_PER_SEC;
li_depois =  cplex_denso.getObjValue();
melhor_LB = cplex_denso.getObjValue();
fim_CPU = clock();
cout << "New LB: " << li_depois << endl;
facilidades_cheias.clear();
facilidades_fracionarias.clear();
facilidades_fechadas.clear();
soma_capacidades_cheias = 0;
  for (int j = 0; j < m; j++)
  {
      
      y_copia[j] = cplex_denso.getValue(y_denso[j])  ;
      if (abs(y_copia[j] - round(y_copia[j]) ) < 0.00001) 
      {
      y_copia[j] = round(y_copia[j]);
          facilidades_cheias.push_back(j);
             soma_capacidades_cheias = soma_capacidades_cheias + p[j];
      }
      else
      {
            facilidades_fracionarias.push_back(j);
      }

  }
if (facilidades_fracionarias.size() > 0)
{
    solucao_inteira = false;
}
for (int i = 0; i < n; i++)
{
  for (int j = 0; j < m; j++)
  {
      
      x_copia[i][j] = cplex_denso.getValue(x_denso[i][j])  ;
      
      if (  abs(round(x_copia[i][j]) - x_copia[i][j]) >= 0.00001  )
      {
          solucao_inteira = false;
      }
      else
      {
          x_copia[i][j] = round(x_copia[i][j]);
      }
        

  }
  
}
LB = cplex_denso.getObjValue();
    if (solucao_inteira == true)
{
melhor_LB = LB;
cout << "LB - Cutting-Plane: " << setprecision(10) << LB << endl; 
fim_CPU = clock();
double tempo_planos_cortantes = (double) (fim_CPU - inicio_CPU)/CLOCKS_PER_SEC;   
cout << "Integer Solution " << endl;
solucao_CPPIS_inteira = true;
double otimo = LB;
//cout << "Time : "<< (double) (fim_CPU - inicio_CPU)/CLOCKS_PER_SEC << endl;
if (conta_iter == 0)
{
fim_CPU = clock();
tempo_planos_cortantes = (double) (fim_CPU - inicio_CPU_pc_eip)/CLOCKS_PER_SEC;
}
double tempo_planos_cortantes_mais_parcial = (double) (fim_CPU - inicio_CPU_pc_eip)/CLOCKS_PER_SEC;
return;
}
}
else
{
pare_planos_cortantes = true;
}

if (   ( abs((li_depois - li_antes)/(li_antes)) < tolerancia_pc) /* && (conta_implied_bounds + conta_lci + conta_lpi + conta_flow_cover_simples + conta_fenchel_cut < 10*/   )
{
    pare_planos_cortantes = true;
}  
fim_CPU = clock();
if ((double) (fim_CPU - inicio_CPU_pc_eip)/CLOCKS_PER_SEC >= tempo_maximo_pc)
{
    pare_planos_cortantes = true;
}
} /// FIM PARE PLANOS CORTANTES. 
cout << "LB - Cutting-Plane:" << setprecision(10) << cplex_denso.getObjValue()  << endl;
fim_CPU = clock();
double tempo_planos_cortantes = (double) (fim_CPU - inicio_CPU_pc_eip)/CLOCKS_PER_SEC;
cout << endl;
}