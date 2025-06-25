
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
#include <ilcplex/ilocplex.h>
#include "Functions.h"
#include <ctime>

using namespace std;


ILOSTLBEGIN 

void SolveSparseProblem(bool problema_esparso_e_gap, bool solucao_heuristica_e_factivel, int quantidade_iteracoes , vector<int> vetor_locais_promissores, vector<double> f, vector<int> vetor_conjunto_clientes, vector<vector<double>> c, vector<double> p, vector<double> pen, vector<double> d, double &UB, bool coeficientes_inteiros, int tempo_limite_problema_esparso, int &iteracao_que_encontrou_otimo, bool &solucao_heuristica, vector<int> atende, list<vector<double>> &coeficientes_cortes_x, list<vector<int>> &clientes_cortes_x, vector<double> &limitantes_cortes_x, vector<int> &variavel_y_associada_ao_corte, list<vector<double>> &coeficientes_cortes_y, list<vector<int>> &facilidades_cortes_y, vector<double> &limitantes_cortes_y, vector<int> &melhor_atende, int n, int m, int tamanho_minimo, double &lb_esparso, double &LB_modelo_otimo)
{
// primeiro_UB_factivel
cout << "Solve Sparse Problem" << endl;
IloEnv env;
IloModel mod(env);
IloCplex cplex(mod);

if ( (problema_esparso_e_gap == true) && (quantidade_iteracoes == 1 ) && ( vetor_locais_promissores.size() <= tamanho_minimo ) )
{

double soma_custos_fixos = 0;
for (int j = 0; j < vetor_locais_promissores.size(); j++)
{
   soma_custos_fixos = soma_custos_fixos + f[vetor_locais_promissores[j]];
}
;
IloArray<IloNumVarArray> x1(env, vetor_conjunto_clientes.size());
for(int i = 0; i < vetor_conjunto_clientes.size(); i++) {
    x1[i] = IloNumVarArray(env, vetor_locais_promissores.size(), 0, 1, ILOBOOL);
}

IloNumVarArray h1(env, vetor_locais_promissores.size(), 0, IloInfinity, ILOFLOAT);
IloExpr expfo(env);
for (int j = 0; j <  vetor_locais_promissores.size(); j++){
                expfo += pow(10,5)*h1[j];
		        for(int i = 0; i < vetor_conjunto_clientes.size(); i++){
		            expfo += c[vetor_conjunto_clientes[i]][vetor_locais_promissores[j]]* d[vetor_conjunto_clientes[i]] * x1[i][j];
            	}
                }
            expfo = expfo + soma_custos_fixos;
            IloAdd(mod, IloMinimize(env, expfo));
expfo.end();

            


for (int j = 0; j < vetor_locais_promissores.size(); j++){
			IloExpr r3(env);
			for (int i = 0; i < vetor_conjunto_clientes.size(); i++){
				r3 += d[vetor_conjunto_clientes[i]]*x1[i][j];
			}
            r3 = r3 - h1[j];
			mod.add(r3  <= p[vetor_locais_promissores[j]] );  //tanto que cabe em cada uma
			r3.end();    
}


//restricao de indivisibilidade*****************************************************************
for (int i = 0; i < vetor_conjunto_clientes.size(); i++){
			IloExpr r2(env);
			for (int j = 0; j < vetor_locais_promissores.size(); j++){
				r2 += x1[i][j];
			}
         	mod.add(r2 == 1);
			r2.end();
    }

            
// Não adicioana corte fo <= UB          

// primeiro UB factivel é false: É necessário adicionar solucao inicial para o solver


/*
vector<vector<int>> start1(n, vector<int>(vetor_locais_promissores.size()));
for (int i = 0; i < n; i++)
{
      for (int j = 0; j < vetor_locais_promissores.size(); j++)
      {
          start1[i][j] = 0;
          if (atende[i] == vetor_locais_promissores[j])
          {
              start1[i][j] = 1;
              //aux_y[vetor_locais_promissores[j]] = 1;
              //fo_aloc = fo_aloc + c[i][vetor_locais_promissores[j]]*d[i];

          }
      }
}
IloNumVarArray startVar(env);
IloNumArray startVal(env);
for (int i = 0; i < n; i++){
         for (int j = 0; j < vetor_locais_promissores.size(); j++) {
             startVar.add(x1[i][j]);
             startVal.add(start1[i][j]);
         }
}
     cplex.addMIPStart(startVar, startVal);
     startVal.end();
     startVar.end();

*/


//=================== configuracoes do solver =============================== 

///////////////////////////////////////////////////////////////////////
////////////ADICIONA CORTES PRODUZIDOS NO NÓ RAIZ SEM VARIÁVEL Y //////

int quantidade_cortes = clientes_cortes_x.size();
int iterador_limitante = 0;
auto it1 = clientes_cortes_x.begin();
auto it2 = coeficientes_cortes_x.begin();


for (int j = 0; j < variavel_y_associada_ao_corte.size(); j++)
{
   /// qual é a posicao que variavel_y_associada_ao_corte_ocupa_em_facilidades_proximas;
   int posicao = -1;
   for (int k = 0; k < vetor_locais_promissores.size(); k++)
   {
      if (variavel_y_associada_ao_corte[j] == vetor_locais_promissores[k])
      {
          posicao = k;
          break;
      }
   }
   if (posicao < 0) 
   {
         std::advance(it1, 1);
         std::advance(it2, 1);
         iterador_limitante = iterador_limitante + 1; 
         continue; 
   }
   


    vector<int> clientes_corte;
    for (int elemento : *it1) {
     clientes_corte.push_back(elemento);
    }
    vector<double> coeficientes_corte;
    for (double elemento : *it2) {
    coeficientes_corte.push_back(elemento);
    }
    IloExpr r2(env);
    for (int i = 0; i < clientes_corte.size(); i++)
    {
       r2 += x1[clientes_corte[i]][posicao]*coeficientes_corte[i];
    }
    mod.add(r2 <= limitantes_cortes_x[iterador_limitante]);
	r2.end();

    std::advance(it1, 1);
     std::advance(it2, 1);
     iterador_limitante = iterador_limitante + 1; 
}
    
////////////////////////////////passa solucao corrente //////////////////////


cplex.setWarning(env.getNullStream()); // Eliminar warnings
//cplex.setOut(env.getNullStream()); /// Eliminar os logs do solver
cplex.setParam(IloCplex::ParallelMode, 1);
cplex.setParam(IloCplex::Threads,1);
cplex.setParam(IloCplex::Param::MIP::Strategy::File,3);

if (coeficientes_inteiros == true)
{
cplex.setParam(IloCplex::Param::MIP::Tolerances::AbsMIPGap, 0.99);
cplex.setParam(IloCplex::EpGap, 0);
}

cplex.setParam(IloCplex::TiLim,tempo_limite_problema_esparso);
cplex.solve();
cout << "status cplex: " <<  cplex.getStatus() << endl;

if  (cplex.getStatus() == IloAlgorithm::Optimal)  
{
    if ( cplex.getObjValue() < UB)
    {
    UB = cplex.getObjValue(); /// atualiza UB
    LB_modelo_otimo = cplex.getBestObjValue();   

    cout << "Uptade UB : " << UB << endl;

    iteracao_que_encontrou_otimo = quantidade_iteracoes;

    vector<double> penalidade_auxiliar(m);
    for (int j = 0; j < m; j++)
    {
       penalidade_auxiliar[j] = -p[j];
    }
    

    for (int i = 0; i < n; i++)
    {

       for (int j = 0; j < vetor_locais_promissores.size(); j++)
       {
       if (cplex.getValue(x1[i][j]) > 0.9)
        {
          melhor_atende[i] = vetor_locais_promissores[j];
          penalidade_auxiliar[vetor_locais_promissores[j]] = penalidade_auxiliar[vetor_locais_promissores[j]] + d[i];
          break;
        }


        }

    }
    /*
    iter_maior_que_um = true;
    for (int j = 0; j < m; j++)
    {
        if (penalidade_auxiliar[j] > 0)
        {
            iter_maior_que_um = false;
            break;
        }
    }
    */

    }
    if ( cplex.getBestObjValue() < lb_esparso)
    {
        lb_esparso = cplex.getBestObjValue();
    }


}
if   (cplex.getStatus() == IloAlgorithm::Unknown) 
{
  solucao_heuristica = true;

}
if (cplex.getStatus() == IloAlgorithm::Feasible)
{

    if ( cplex.getObjValue() < UB)
    {
    UB = cplex.getObjValue();
   
        LB_modelo_otimo = cplex.getBestObjValue();  
    cout << "Uptade UB : " << UB << endl;

    iteracao_que_encontrou_otimo = quantidade_iteracoes;

    vector<double> penalidade_auxiliar(m);
    for (int j = 0; j < m; j++)
    {
       penalidade_auxiliar[j] = -p[j];
    }
    

    for (int i = 0; i < n; i++)
    {

       for (int j = 0; j < vetor_locais_promissores.size(); j++)
       {
       if (cplex.getValue(x1[i][j]) > 0.9)
        {
          melhor_atende[i] = vetor_locais_promissores[j];
          penalidade_auxiliar[vetor_locais_promissores[j]] = penalidade_auxiliar[vetor_locais_promissores[j]] + d[i];
          break;
        }


        }

    }
    /*
    iter_maior_que_um = true;
    for (int j = 0; j < m; j++)
    {
        if (penalidade_auxiliar[j] > 0)
        {
            iter_maior_que_um = false;
            break;
        }
    }
    */

    }
        if ( cplex.getBestObjValue() < lb_esparso)
    {
        lb_esparso = cplex.getBestObjValue();
    }

solucao_heuristica = true;

}




return;
} 
if ( (problema_esparso_e_gap == false) && (quantidade_iteracoes == 1) && (vetor_locais_promissores.size() <= tamanho_minimo) )
{



//IloEnv env;
//IloModel mod(env);
//IloCplex cplex(mod);
IloNumVarArray y(env, vetor_locais_promissores.size(), 0, 1, ILOBOOL);
IloArray<IloNumVarArray> x1(env, vetor_conjunto_clientes.size());
//IloArray<IloNumVarArray> x2(env, vetor_conjunto_clientes.size());

for(int i = 0; i < vetor_conjunto_clientes.size(); i++) {
    x1[i] = IloNumVarArray(env, vetor_locais_promissores.size(), 0, 1, ILOBOOL);
}


IloExpr expfo(env);
for (int j = 0; j <  vetor_locais_promissores.size(); j++){
	expfo += f[vetor_locais_promissores[j]] * y[j];
}
for (int j = 0; j <  vetor_locais_promissores.size(); j++){
		for(int i = 0; i < vetor_conjunto_clientes.size(); i++){
		expfo += c[vetor_conjunto_clientes[i]][vetor_locais_promissores[j]]* d[vetor_conjunto_clientes[i]] * x1[i][j];
	}
}

IloNumVarArray h1(env, vetor_locais_promissores.size(), 0, IloInfinity, ILOFLOAT);

for (int j = 0; j <  vetor_locais_promissores.size(); j++){
	expfo += 100000*h1[j];
}

IloAdd(mod, IloMinimize(env, expfo));
expfo.end();
            


			for (int j = 0; j < vetor_locais_promissores.size(); j++){
			IloExpr r3(env);
			for (int i = 0; i < vetor_conjunto_clientes.size(); i++){
				r3 += d[vetor_conjunto_clientes[i]]*x1[i][j];
			}
            r3 = r3 - h1[j];
			mod.add(r3  <= p[vetor_locais_promissores[j]]*y[j]);  //tanto que cabe em cada uma
			r3.end();
		    }
            

//restricao de indivisibilidade*****************************************************************
for (int i = 0; i < vetor_conjunto_clientes.size(); i++){
			IloExpr r2(env);
			for (int j = 0; j < vetor_locais_promissores.size(); j++){
				r2 += x1[i][j];
			}
         	mod.add(r2 == 1);
			r2.end();
    }

            
// Não adicioana corte fo <= UB          

// primeiro UB factivel é false: É necessário adicionar solucao inicial para o solver




///////////////////////////////////////////////////////////////////////
////////////ADICIONA CORTES PRODUZIDOS NO NÓ RAIZ COM VARIÁVEL Y //////

int quantidade_cortes = clientes_cortes_x.size();
int iterador_limitante = 0;
auto it1 = clientes_cortes_x.begin();
auto it2 = coeficientes_cortes_x.begin();


for (int j = 0; j < variavel_y_associada_ao_corte.size(); j++)
{
   /// qual é a posicao que variavel_y_associada_ao_corte_ocupa_em_facilidades_proximas;
   int posicao = -1;
   for (int k = 0; k < vetor_locais_promissores.size(); k++)
   {
      if (variavel_y_associada_ao_corte[j] == vetor_locais_promissores[k])
      {
          posicao = k;
          break;
      }
   }
   if (posicao < 0) 
   {
         std::advance(it1, 1);
         std::advance(it2, 1);
         iterador_limitante = iterador_limitante + 1; 
         continue; 
   }
   


    vector<int> clientes_corte;
    for (int elemento : *it1) {
     clientes_corte.push_back(elemento);
    }
    vector<double> coeficientes_corte;
    for (double elemento : *it2) {
    coeficientes_corte.push_back(elemento);
    }
    IloExpr r2(env);
    for (int i = 0; i < clientes_corte.size(); i++)
    {
       r2 += x1[clientes_corte[i]][posicao]*coeficientes_corte[i];
    }
    mod.add(r2 <= limitantes_cortes_x[iterador_limitante]*y[posicao]);
	r2.end();

    std::advance(it1, 1);
     std::advance(it2, 1);
     iterador_limitante = iterador_limitante + 1; 
}
    
////////////////////////////////passa solucao corrente //////////////////////


//=================== configuracoes do solver =============================== 

cplex.setWarning(env.getNullStream()); // Eliminar warnings
//cplex.setOut(env.getNullStream()); /// Eliminar os logs do solver
cplex.setParam(IloCplex::ParallelMode, 1);
cplex.setParam(IloCplex::Threads,1);
cplex.setParam(IloCplex::Param::MIP::Strategy::File,3);

if (coeficientes_inteiros == true)
{
cplex.setParam(IloCplex::Param::MIP::Tolerances::AbsMIPGap, 0.99);
cplex.setParam(IloCplex::EpGap,0);
}

cplex.setParam(IloCplex::TiLim,tempo_limite_problema_esparso);
cplex.solve();
cout << "status cplex: " <<  cplex.getStatus() << endl;

if  (cplex.getStatus() == IloAlgorithm::Optimal)  
{
    if ( cplex.getObjValue() < UB)
    {
    UB = cplex.getObjValue();
        LB_modelo_otimo = cplex.getBestObjValue();  
    cout << "Uptade UB : " << UB << endl;

    iteracao_que_encontrou_otimo = quantidade_iteracoes;

    vector<double> penalidade_auxiliar(m);
    for (int j = 0; j < m; j++)
    {
       penalidade_auxiliar[j] = -p[j];
    }
    

    for (int i = 0; i < n; i++)
    {

       for (int j = 0; j < vetor_locais_promissores.size(); j++)
       {
       if (cplex.getValue(x1[i][j]) > 0.9)
        {
          melhor_atende[i] = vetor_locais_promissores[j];
          penalidade_auxiliar[vetor_locais_promissores[j]] = penalidade_auxiliar[vetor_locais_promissores[j]] + d[i];
          break;
        }


        }

    }
    /*
    iter_maior_que_um = true;
    for (int j = 0; j < m; j++)
    {
        if (penalidade_auxiliar[j] > 0)
        {
            iter_maior_que_um = false;
            break;
        }
    }
    */

    }
        if ( cplex.getBestObjValue() < lb_esparso)
    {
        lb_esparso = cplex.getBestObjValue();
    }



}
if   (cplex.getStatus() == IloAlgorithm::Unknown) 
{
  solucao_heuristica = true;

}
if (cplex.getStatus() == IloAlgorithm::Feasible)
{


    if ( cplex.getObjValue() < UB)
    {
    UB = cplex.getObjValue();
        LB_modelo_otimo = cplex.getBestObjValue();  
    cout << "Uptade UB : " << UB << endl;

    iteracao_que_encontrou_otimo = quantidade_iteracoes;

    vector<double> penalidade_auxiliar(m);
    for (int j = 0; j < m; j++)
    {
       penalidade_auxiliar[j] = -p[j];
    }
    

    for (int i = 0; i < n; i++)
    {

       for (int j = 0; j < vetor_locais_promissores.size(); j++)
       {
       if (cplex.getValue(x1[i][j]) > 0.9)
        {
          melhor_atende[i] = vetor_locais_promissores[j];
          penalidade_auxiliar[vetor_locais_promissores[j]] = penalidade_auxiliar[vetor_locais_promissores[j]] + d[i];
          break;
        }


        }

    }
    /*
    iter_maior_que_um = true;
    for (int j = 0; j < m; j++)
    {
        if (penalidade_auxiliar[j] > 0)
        {
            iter_maior_que_um = false;
            break;
        }
    }
    */

    }
        if ( cplex.getBestObjValue() < lb_esparso)
    {
        lb_esparso = cplex.getBestObjValue();
    }

solucao_heuristica = true;

}


return;

}

if ( (problema_esparso_e_gap == true) && (quantidade_iteracoes > 1) )
{



double soma_custos_fixos = 0;
for (int j = 0; j < vetor_locais_promissores.size(); j++)
{
   soma_custos_fixos = soma_custos_fixos + f[vetor_locais_promissores[j]];
}
//IloEnv env;
//IloModel mod(env);
//IloCplex cplex(mod);
IloArray<IloNumVarArray> x1(env, vetor_conjunto_clientes.size());
for(int i = 0; i < vetor_conjunto_clientes.size(); i++) {
    x1[i] = IloNumVarArray(env, vetor_locais_promissores.size(), 0, 1, ILOBOOL);
}

// IloNumVarArray h1(env, vetor_locais_promissores.size(), 0, IloInfinity, ILOFLOAT);
IloExpr expfo(env);
for (int j = 0; j <  vetor_locais_promissores.size(); j++){
		        for(int i = 0; i < vetor_conjunto_clientes.size(); i++){
		            expfo += c[vetor_conjunto_clientes[i]][vetor_locais_promissores[j]]* d[vetor_conjunto_clientes[i]] * x1[i][j];
            	}
                }
            expfo = expfo + soma_custos_fixos;
            IloAdd(mod, IloMinimize(env, expfo));
expfo.end();

            


for (int j = 0; j < vetor_locais_promissores.size(); j++){
			IloExpr r3(env);
			for (int i = 0; i < vetor_conjunto_clientes.size(); i++){
				r3 += d[vetor_conjunto_clientes[i]]*x1[i][j];
			}
            //r3 = r3 - h1[j];
			mod.add(r3  <= p[vetor_locais_promissores[j]] );  //tanto que cabe em cada uma
			r3.end();    
}


//restricao de indivisibilidade*****************************************************************
for (int i = 0; i < vetor_conjunto_clientes.size(); i++){
			IloExpr r2(env);
			for (int j = 0; j < vetor_locais_promissores.size(); j++){
				r2 += x1[i][j];
			}
         	mod.add(r2 == 1);
			r2.end();
    }

            
//adicioana corte fo <= UB          
/*
IloExpr r15(env);
    for (int j = 0; j <  vetor_locais_promissores.size(); j++){
		    for(int i = 0; i < vetor_conjunto_clientes.size(); i++){
		    r15 += c[vetor_conjunto_clientes[i]][vetor_locais_promissores[j]]* d[vetor_conjunto_clientes[i]] * x1[i][j];
	        }
            }
        r15 = r15 + soma_custos_fixos;
        //cout << "valor de UB:" << UB << endl;
        if (coeficientes_inteiros == true)
            {   
                mod.add(r15 <= UB - 1);
            }
        else
            {
                mod.add(r15 <= UB - 0.001);    
            }
            */

// primeiro UB factivel: não é necessário adicionar solucao inicial para o solver

//=================== configuracoes do solver =============================== 




///////////////////////////////////////////////////////////////////////
////////////ADICIONA CORTES PRODUZIDOS NO NÓ RAIZ SEM VARIÁVEL Y //////

int quantidade_cortes = clientes_cortes_x.size();
int iterador_limitante = 0;
auto it1 = clientes_cortes_x.begin();
auto it2 = coeficientes_cortes_x.begin();


for (int j = 0; j < variavel_y_associada_ao_corte.size(); j++)
{
   /// qual é a posicao que variavel_y_associada_ao_corte_ocupa_em_facilidades_proximas;
   int posicao = -1;
   for (int k = 0; k < vetor_locais_promissores.size(); k++)
   {
      if (variavel_y_associada_ao_corte[j] == vetor_locais_promissores[k])
      {
          posicao = k;
          break;
      }
   }
   if (posicao < 0) 
   {
         std::advance(it1, 1);
         std::advance(it2, 1);
         iterador_limitante = iterador_limitante + 1; 
         continue; 
   }
   


    vector<int> clientes_corte;
    for (int elemento : *it1) {
     clientes_corte.push_back(elemento);
    }
    vector<double> coeficientes_corte;
    for (double elemento : *it2) {
    coeficientes_corte.push_back(elemento);
    }
    IloExpr r2(env);
    for (int i = 0; i < clientes_corte.size(); i++)
    {
       r2 += x1[clientes_corte[i]][posicao]*coeficientes_corte[i];
    }
    mod.add(r2 <= limitantes_cortes_x[iterador_limitante]);
	r2.end();

    std::advance(it1, 1);
     std::advance(it2, 1);
     iterador_limitante = iterador_limitante + 1; 
}
    
////////////////////////////////passa solucao corrente //////////////////////


cplex.setWarning(env.getNullStream()); // Eliminar warnings
//cplex.setOut(env.getNullStream()); /// Eliminar os logs do solver
cplex.setParam(IloCplex::ParallelMode, 1);
cplex.setParam(IloCplex::Threads,1);
cplex.setParam(IloCplex::Param::MIP::Strategy::File,3);

if (coeficientes_inteiros == true)
{
cplex.setParam(IloCplex::CutUp, UB - 1);    
cplex.setParam(IloCplex::Param::MIP::Tolerances::AbsMIPGap, 0.99);
cplex.setParam(IloCplex::EpGap,0);
}
else
{
    cplex.setParam(IloCplex::CutUp, UB - 0.001); 
}

cplex.setParam(IloCplex::TiLim,tempo_limite_problema_esparso);
cplex.solve();
cout << "status cplex: " <<  cplex.getStatus() << endl;
if  (cplex.getStatus() == IloAlgorithm::Optimal)  // encontrou solucao melhor do que UB corrente. mas UB corrente é factivel
{

    UB = cplex.getObjValue();
        LB_modelo_otimo = cplex.getBestObjValue();  
    cout << "Uptade UB : " << UB << endl;
    //primeira_solucao_factivel_encontrada = factivel;
    //cout << "novo apos resolucao do modelo: " << UB << endl;
    iteracao_que_encontrou_otimo = quantidade_iteracoes;
    //cout << "limitante inferior do último problema: " << cplex.getBestObjValue() << endl;
    //arq_saida1 <<  " Novo UB : " << fo << " tempo cplex para otimilidade subproblema esparso :" << (double) (fim_CPU5 - inicio_CPU_cut_and_solve5)/CLOCKS_PER_SEC << endl;
    /*
    for (int j = 0; j < m; j++)
    {
       pen_melhor_UB[j] = -p[j];
    }
    */
    for (int i = 0; i < n; i++)
    {

       for (int j = 0; j < vetor_locais_promissores.size(); j++)
       {
       if (cplex.getValue(x1[i][j]) > 0.9)
        {
          melhor_atende[i] = vetor_locais_promissores[j];
          //pen_melhor_UB[vetor_locais_promissores[j]] = pen_melhor_UB[vetor_locais_promissores[j]] + d[i];
          break;
        }


        }

    }


        if ( cplex.getBestObjValue() < lb_esparso)
    {
        lb_esparso = cplex.getBestObjValue();
    }


    


}
if   (cplex.getStatus() == IloAlgorithm::Unknown) 
{
  solucao_heuristica = true;

}
if (cplex.getStatus() == IloAlgorithm::Feasible)
{

if (cplex.getObjValue() < UB)
{
    UB = cplex.getObjValue();
        LB_modelo_otimo = cplex.getBestObjValue();  
    //atualiza melhor atende
        for (int i = 0; i < n; i++)
    {

       for (int j = 0; j < vetor_locais_promissores.size(); j++)
       {
       if (cplex.getValue(x1[i][j]) > 0.9)
        {
          melhor_atende[i] = vetor_locais_promissores[j];
          break;
        }


        }

    }
iteracao_que_encontrou_otimo = quantidade_iteracoes;
}
solucao_heuristica = true;
    if ( cplex.getBestObjValue() < lb_esparso)
    {
        lb_esparso = cplex.getBestObjValue();
    }

}



return;
}


if ( (problema_esparso_e_gap == true) && (quantidade_iteracoes == 1 )  )
{
// como primeiro ub é factivel entao rodou a heuristica 
// solucao_heuristica_e_factivel == true: não precisa criar variáveis de excesso 
double soma_custos_fixos = 0;
for (int j = 0; j < vetor_locais_promissores.size(); j++)
{
   soma_custos_fixos = soma_custos_fixos + f[vetor_locais_promissores[j]];
}
//IloEnv env;
//IloModel mod(env);
//IloCplex cplex(mod);
IloArray<IloNumVarArray> x1(env, vetor_conjunto_clientes.size());
for(int i = 0; i < vetor_conjunto_clientes.size(); i++) {
    x1[i] = IloNumVarArray(env, vetor_locais_promissores.size(), 0, 1, ILOBOOL);
}

// IloNumVarArray h1(env, vetor_locais_promissores.size(), 0, IloInfinity, ILOFLOAT);
IloExpr expfo(env);
for (int j = 0; j <  vetor_locais_promissores.size(); j++){
		        for(int i = 0; i < vetor_conjunto_clientes.size(); i++){
		            expfo += c[vetor_conjunto_clientes[i]][vetor_locais_promissores[j]]* d[vetor_conjunto_clientes[i]] * x1[i][j];
            	}
                }
            expfo = expfo + soma_custos_fixos;
            IloAdd(mod, IloMinimize(env, expfo));
expfo.end();

            


for (int j = 0; j < vetor_locais_promissores.size(); j++){
			IloExpr r3(env);
			for (int i = 0; i < vetor_conjunto_clientes.size(); i++){
				r3 += d[vetor_conjunto_clientes[i]]*x1[i][j];
			}
            //r3 = r3 - h1[j];
			mod.add(r3  <= p[vetor_locais_promissores[j]] );  //tanto que cabe em cada uma
			r3.end();    
}


//restricao de indivisibilidade*****************************************************************
for (int i = 0; i < vetor_conjunto_clientes.size(); i++){
			IloExpr r2(env);
			for (int j = 0; j < vetor_locais_promissores.size(); j++){
				r2 += x1[i][j];
			}
         	mod.add(r2 == 1);
			r2.end();
    }

            
// Não adicioana corte fo <= UB          

// primeiro UB factivel é false: É necessário adicionar solucao inicial para o solver



/*
vector<vector<int>> start1(n, vector<int>(vetor_locais_promissores.size()));
for (int i = 0; i < n; i++)
{
      for (int j = 0; j < vetor_locais_promissores.size(); j++)
      {
          start1[i][j] = 0;
          if (atende[i] == vetor_locais_promissores[j])
          {
              start1[i][j] = 1;
              //aux_y[vetor_locais_promissores[j]] = 1;
              //fo_aloc = fo_aloc + c[i][vetor_locais_promissores[j]]*d[i];

          }
      }
}
IloNumVarArray startVar(env);
IloNumArray startVal(env);
for (int i = 0; i < n; i++){
         for (int j = 0; j < vetor_locais_promissores.size(); j++) {
             startVar.add(x1[i][j]);
             startVal.add(start1[i][j]);
         }
}
     cplex.addMIPStart(startVar, startVal);
     startVal.end();
     startVar.end();
*/

//=================== configuracoes do solver =============================== 

///////////////////////////////////////////////////////////////////////
////////////ADICIONA CORTES PRODUZIDOS NO NÓ RAIZ SEM VARIÁVEL Y //////



int quantidade_cortes = clientes_cortes_x.size();
int iterador_limitante = 0;
auto it1 = clientes_cortes_x.begin();
auto it2 = coeficientes_cortes_x.begin();


for (int j = 0; j < variavel_y_associada_ao_corte.size(); j++)
{
   /// qual é a posicao que variavel_y_associada_ao_corte_ocupa_em_facilidades_proximas;
   int posicao = -1;
   for (int k = 0; k < vetor_locais_promissores.size(); k++)
   {
      if (variavel_y_associada_ao_corte[j] == vetor_locais_promissores[k])
      {
          posicao = k;
          break;
      }
   }
   if (posicao < 0) 
   {
         std::advance(it1, 1);
         std::advance(it2, 1);
         iterador_limitante = iterador_limitante + 1; 
         continue; 
   }
   


    vector<int> clientes_corte;
    for (int elemento : *it1) {
     clientes_corte.push_back(elemento);
    }
    vector<double> coeficientes_corte;
    for (double elemento : *it2) {
    coeficientes_corte.push_back(elemento);
    }
    IloExpr r2(env);
    for (int i = 0; i < clientes_corte.size(); i++)
    {
       r2 += x1[clientes_corte[i]][posicao]*coeficientes_corte[i];
    }
    mod.add(r2 <= limitantes_cortes_x[iterador_limitante]);
	r2.end();

    std::advance(it1, 1);
     std::advance(it2, 1);
     iterador_limitante = iterador_limitante + 1; 
}

   

////////////////////////////////passa solucao corrente //////////////////////


/*
IloExpr r15(env);
    for (int j = 0; j <  vetor_locais_promissores.size(); j++){
		    for(int i = 0; i < vetor_conjunto_clientes.size(); i++){
		    r15 += c[vetor_conjunto_clientes[i]][vetor_locais_promissores[j]]* d[vetor_conjunto_clientes[i]] * x1[i][j];
	        }
            }
        r15 = r15 + soma_custos_fixos;
        //cout << "valor de UB:" << UB << endl;
        if (coeficientes_inteiros == true)
            {   
                mod.add(r15 <= UB - 1);
            }
        else
            {
                mod.add(r15 <= UB - 0.001);    
            }
*/


cplex.setWarning(env.getNullStream()); // Eliminar warnings
//cplex.setOut(env.getNullStream()); /// Eliminar os logs do solver
cplex.setParam(IloCplex::ParallelMode, 1);
cplex.setParam(IloCplex::Threads,1);
cplex.setParam(IloCplex::Param::MIP::Strategy::File,3);

if (coeficientes_inteiros == true)
{
cplex.setParam(IloCplex::CutUp, UB - 1);    
cplex.setParam(IloCplex::Param::MIP::Tolerances::AbsMIPGap, 0.99);
cplex.setParam(IloCplex::EpGap,0);
}
else
{
    cplex.setParam(IloCplex::CutUp, UB - 0.001); 
}



cplex.setParam(IloCplex::TiLim,tempo_limite_problema_esparso);
cplex.solve();
cout << "status cplex: " <<  cplex.getStatus() << endl;
if  (cplex.getStatus() == IloAlgorithm::Optimal)  
{
    if ( cplex.getObjValue() < UB)
    {
    UB = cplex.getObjValue();
        LB_modelo_otimo = cplex.getBestObjValue();  
    cout << "Uptade UB : " << UB << endl;

    iteracao_que_encontrou_otimo = quantidade_iteracoes;

    vector<double> penalidade_auxiliar(m);
    for (int j = 0; j < m; j++)
    {
       penalidade_auxiliar[j] = -p[j];
    }
    

    for (int i = 0; i < n; i++)
    {

       for (int j = 0; j < vetor_locais_promissores.size(); j++)
       {
       if (cplex.getValue(x1[i][j]) > 0.9)
        {
          melhor_atende[i] = vetor_locais_promissores[j];
          penalidade_auxiliar[vetor_locais_promissores[j]] = penalidade_auxiliar[vetor_locais_promissores[j]] + d[i];
          break;
        }


        }

    }
    /*
    iter_maior_que_um = true;
    for (int j = 0; j < m; j++)
    {
        if (penalidade_auxiliar[j] > 0)
        {
            iter_maior_que_um = false;
            break;
        }
    }
    */

    }
    
    if ( cplex.getBestObjValue() < lb_esparso)
    {
        lb_esparso = cplex.getBestObjValue();
    }


}
if   (cplex.getStatus() == IloAlgorithm::Unknown) 
{
  solucao_heuristica = true;

}
if (cplex.getStatus() == IloAlgorithm::Feasible)
{

    if ( cplex.getObjValue() < UB)
    {
    UB = cplex.getObjValue();
        LB_modelo_otimo = cplex.getBestObjValue();  
    cout << "Uptade UB : " << UB << endl;

    iteracao_que_encontrou_otimo = quantidade_iteracoes;

    vector<double> penalidade_auxiliar(m);
    for (int j = 0; j < m; j++)
    {
       penalidade_auxiliar[j] = -p[j];
    }
    

    for (int i = 0; i < n; i++)
    {

       for (int j = 0; j < vetor_locais_promissores.size(); j++)
       {
       if (cplex.getValue(x1[i][j]) > 0.9)
        {
          melhor_atende[i] = vetor_locais_promissores[j];
          penalidade_auxiliar[vetor_locais_promissores[j]] = penalidade_auxiliar[vetor_locais_promissores[j]] + d[i];
          break;
        }


        }

    }
    /*
    iter_maior_que_um = true;
    for (int j = 0; j < m; j++)
    {
        if (penalidade_auxiliar[j] > 0)
        {
            iter_maior_que_um = false;
            break;
        }
    }
    */

    }
solucao_heuristica = true;
    if ( cplex.getBestObjValue() < lb_esparso)
    {
        lb_esparso = cplex.getBestObjValue();
    }

}


return;

}


if ( (problema_esparso_e_gap == false) && (quantidade_iteracoes > 1)   )
{

IloNumVarArray y(env, vetor_locais_promissores.size(), 0, 1, ILOBOOL);
IloArray<IloNumVarArray> x1(env, vetor_conjunto_clientes.size());

for(int i = 0; i < vetor_conjunto_clientes.size(); i++) {
    x1[i] = IloNumVarArray(env, vetor_locais_promissores.size(), 0, 1, ILOBOOL);
}
IloExpr expfo(env);
for (int j = 0; j <  vetor_locais_promissores.size(); j++){
	expfo += f[vetor_locais_promissores[j]] * y[j];
}
for (int j = 0; j <  vetor_locais_promissores.size(); j++){
		for(int i = 0; i < vetor_conjunto_clientes.size(); i++){
		expfo += c[vetor_conjunto_clientes[i]][vetor_locais_promissores[j]]* d[vetor_conjunto_clientes[i]] * x1[i][j];
	}
}

IloAdd(mod, IloMinimize(env, expfo));
expfo.end();
            

			for (int j = 0; j < vetor_locais_promissores.size(); j++){
			IloExpr r3(env);
			for (int i = 0; i < vetor_conjunto_clientes.size(); i++){
				r3 += d[vetor_conjunto_clientes[i]]*x1[i][j];
			}
            //r3 = r3 - h1[j];
			mod.add(r3  <= p[vetor_locais_promissores[j]]*y[j]);  //tanto que cabe em cada uma
			r3.end();
		    }
            
//restricao de indivisibilidade*****************************************************************
for (int i = 0; i < vetor_conjunto_clientes.size(); i++){
			IloExpr r2(env);
			for (int j = 0; j < vetor_locais_promissores.size(); j++){
				r2 += x1[i][j];
			}
         	mod.add(r2 == 1);
			r2.end();
    }

////////////ADICIONA CORTES PRODUZIDOS NO NÓ RAIZ COM VARIÁVEL Y //////
int quantidade_cortes = clientes_cortes_x.size();
int iterador_limitante = 0;
auto it1 = clientes_cortes_x.begin();
auto it2 = coeficientes_cortes_x.begin();
for (int j = 0; j < variavel_y_associada_ao_corte.size(); j++)
{
   /// qual é a posicao que variavel_y_associada_ao_corte_ocupa_em_facilidades_proximas;
   int posicao = -1;
   for (int k = 0; k < vetor_locais_promissores.size(); k++)
   {
      if (variavel_y_associada_ao_corte[j] == vetor_locais_promissores[k])
      {
          posicao = k;
          break;
      }
   }
   if (posicao < 0) 
   {
         std::advance(it1, 1);
         std::advance(it2, 1);
         iterador_limitante = iterador_limitante + 1; 
         continue; 
   }
      vector<int> clientes_corte;
    for (int elemento : *it1) {
     clientes_corte.push_back(elemento);
    }
    vector<double> coeficientes_corte;
    for (double elemento : *it2) {
    coeficientes_corte.push_back(elemento);
    }
    IloExpr r2(env);
    for (int i = 0; i < clientes_corte.size(); i++)
    {
       r2 += x1[clientes_corte[i]][posicao]*coeficientes_corte[i];
    }
    mod.add(r2 <= limitantes_cortes_x[iterador_limitante]*y[posicao]);
	r2.end();

    std::advance(it1, 1);
     std::advance(it2, 1);
     iterador_limitante = iterador_limitante + 1; 
}
    
////////////////////////////////passa solucao corrente //////////////////////
cplex.setWarning(env.getNullStream()); // Eliminar warnings
//cplex.setOut(env.getNullStream()); /// Eliminar os logs do solver
cplex.setParam(IloCplex::ParallelMode, 1);
cplex.setParam(IloCplex::Threads,1);
cplex.setParam(IloCplex::Param::MIP::Strategy::File,3);

if (coeficientes_inteiros == true)
{
cplex.setParam(IloCplex::CutUp, UB - 1);    
cplex.setParam(IloCplex::Param::MIP::Tolerances::AbsMIPGap, 0.99);
cplex.setParam(IloCplex::EpGap,0);
}
else
{
    cplex.setParam(IloCplex::CutUp, UB - 0.001); 
}
cplex.setParam(IloCplex::TiLim,tempo_limite_problema_esparso);
cplex.solve();
cout << "status cplex: " <<  cplex.getStatus() << endl;
if  (cplex.getStatus() == IloAlgorithm::Optimal)  // encontrou solucao melhor do que UB corrente. mas UB corrente é factivel
{

    UB = cplex.getObjValue();
        LB_modelo_otimo = cplex.getBestObjValue();  
    cout << "Uptade UB : " << UB << endl;
    iteracao_que_encontrou_otimo = quantidade_iteracoes;
    for (int i = 0; i < n; i++)
    {

       for (int j = 0; j < vetor_locais_promissores.size(); j++)
       {
       if (cplex.getValue(x1[i][j]) > 0.9)
        {
          melhor_atende[i] = vetor_locais_promissores[j];
          break;
        }
        }
    }
    if ( cplex.getBestObjValue() < lb_esparso)
    {
        lb_esparso = cplex.getBestObjValue();
    }
}
if   (cplex.getStatus() == IloAlgorithm::Unknown) 
{
  solucao_heuristica = true;

}
if (cplex.getStatus() == IloAlgorithm::Feasible)
{
if (cplex.getObjValue() < UB)
{
    UB = cplex.getObjValue();
        LB_modelo_otimo = cplex.getBestObjValue();  
        for (int i = 0; i < n; i++)
    {

       for (int j = 0; j < vetor_locais_promissores.size(); j++)
       {
       if (cplex.getValue(x1[i][j]) > 0.9)
        {
          melhor_atende[i] = vetor_locais_promissores[j];
          break;
        }
        }
    }
iteracao_que_encontrou_otimo = quantidade_iteracoes;
}
solucao_heuristica = true;
    if ( cplex.getBestObjValue() < lb_esparso)
    {
        lb_esparso = cplex.getBestObjValue();
    }
}
return;
}
if ( (problema_esparso_e_gap == false) && (quantidade_iteracoes == 1)  )
{
IloNumVarArray y(env, vetor_locais_promissores.size(), 0, 1, ILOBOOL);
IloArray<IloNumVarArray> x1(env, vetor_conjunto_clientes.size());
for(int i = 0; i < vetor_conjunto_clientes.size(); i++) {
    x1[i] = IloNumVarArray(env, vetor_locais_promissores.size(), 0, 1, ILOBOOL);
}
IloExpr expfo(env);
for (int j = 0; j <  vetor_locais_promissores.size(); j++){
	expfo += f[vetor_locais_promissores[j]] * y[j];
}
for (int j = 0; j <  vetor_locais_promissores.size(); j++){
		for(int i = 0; i < vetor_conjunto_clientes.size(); i++){
		expfo += c[vetor_conjunto_clientes[i]][vetor_locais_promissores[j]]* d[vetor_conjunto_clientes[i]] * x1[i][j];
	}
}
IloNumVarArray h1(env, vetor_locais_promissores.size(), 0, IloInfinity, ILOFLOAT);
for (int j = 0; j <  vetor_locais_promissores.size(); j++){
	expfo += 100000*h1[j];
}
IloAdd(mod, IloMinimize(env, expfo));
expfo.end();
			for (int j = 0; j < vetor_locais_promissores.size(); j++){
			IloExpr r3(env);
			for (int i = 0; i < vetor_conjunto_clientes.size(); i++){
				r3 += d[vetor_conjunto_clientes[i]]*x1[i][j];
			}
            r3 = r3 - h1[j];
			mod.add(r3  <= p[vetor_locais_promissores[j]]*y[j]);  //tanto que cabe em cada uma
			r3.end();
		    }
for (int i = 0; i < vetor_conjunto_clientes.size(); i++){
			IloExpr r2(env);
			for (int j = 0; j < vetor_locais_promissores.size(); j++){
				r2 += x1[i][j];
			}
         	mod.add(r2 == 1);
			r2.end();
    }
////////////ADICIONA CORTES PRODUZIDOS NO NÓ RAIZ COM VARIÁVEL Y //////
int quantidade_cortes = clientes_cortes_x.size();
int iterador_limitante = 0;
auto it1 = clientes_cortes_x.begin();
auto it2 = coeficientes_cortes_x.begin();
for (int j = 0; j < variavel_y_associada_ao_corte.size(); j++)
{
   int posicao = -1;
   for (int k = 0; k < vetor_locais_promissores.size(); k++)
   {
      if (variavel_y_associada_ao_corte[j] == vetor_locais_promissores[k])
      {
          posicao = k;
          break;
      }
   }
   if (posicao < 0) 
   {
         std::advance(it1, 1);
         std::advance(it2, 1);
         iterador_limitante = iterador_limitante + 1; 
         continue; 
   }
    vector<int> clientes_corte;
    for (int elemento : *it1) {
     clientes_corte.push_back(elemento);
    }
    vector<double> coeficientes_corte;
    for (double elemento : *it2) {
    coeficientes_corte.push_back(elemento);
    }
    IloExpr r2(env);
    for (int i = 0; i < clientes_corte.size(); i++)
    {
       r2 += x1[clientes_corte[i]][posicao]*coeficientes_corte[i];
    }
    mod.add(r2 <= limitantes_cortes_x[iterador_limitante]*y[posicao]);
	r2.end();
    std::advance(it1, 1);
    std::advance(it2, 1);
     iterador_limitante = iterador_limitante + 1; 
}
cplex.setWarning(env.getNullStream()); // Eliminar warnings
cplex.setParam(IloCplex::ParallelMode, 1);
cplex.setParam(IloCplex::Threads,1);
cplex.setParam(IloCplex::Param::MIP::Strategy::File,3);
if (coeficientes_inteiros == true)
{
cplex.setParam(IloCplex::CutUp, UB - 1);    
cplex.setParam(IloCplex::Param::MIP::Tolerances::AbsMIPGap, 0.99);
cplex.setParam(IloCplex::EpGap,0);
}
else
{
    cplex.setParam(IloCplex::CutUp, UB - 0.001); 
}
cplex.setParam(IloCplex::TiLim,tempo_limite_problema_esparso);
cplex.solve();
cout << "status cplex: " <<  cplex.getStatus() << endl;
if  (cplex.getStatus() == IloAlgorithm::Optimal)  // encontrou solucao melhor do que UB corrente. mas UB corrente é factivel
{
    if (cplex.getObjValue() < UB)
    {
    UB = cplex.getObjValue();
        LB_modelo_otimo = cplex.getBestObjValue();  
    cout << "Uptade UB : " << UB << endl;
    iteracao_que_encontrou_otimo = quantidade_iteracoes;

    for (int i = 0; i < n; i++)
    {
       for (int j = 0; j < vetor_locais_promissores.size(); j++)
       {
       if (cplex.getValue(x1[i][j]) > 0.9)
        {
          melhor_atende[i] = vetor_locais_promissores[j];
          //pen_melhor_UB[vetor_locais_promissores[j]] = pen_melhor_UB[vetor_locais_promissores[j]] + d[i];
          break;
        }
        }
    }
    }
    if ( cplex.getBestObjValue() < lb_esparso)
    {
        lb_esparso = cplex.getBestObjValue();
    }
}
if   (cplex.getStatus() == IloAlgorithm::Unknown) 
{
  solucao_heuristica = true;

}
if (cplex.getStatus() == IloAlgorithm::Feasible)
{

if (cplex.getObjValue() < UB)
{
    UB = cplex.getObjValue();
        LB_modelo_otimo = cplex.getBestObjValue();  
    //atualiza melhor atende
        for (int i = 0; i < n; i++)
    {

       for (int j = 0; j < vetor_locais_promissores.size(); j++)
       {
       if (cplex.getValue(x1[i][j]) > 0.9)
        {
          melhor_atende[i] = vetor_locais_promissores[j];
          break;
        }


        }

    }
iteracao_que_encontrou_otimo = quantidade_iteracoes;
}
solucao_heuristica = true;
    if ( cplex.getBestObjValue() < lb_esparso)
    {
        lb_esparso = cplex.getBestObjValue();
    }

}
return;
}
mod.end();
cplex.end();
env.end();
}

