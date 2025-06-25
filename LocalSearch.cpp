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
#include <ctime>
#include <string.h>
#include <malloc.h>
#include <ilcplex/ilocplex.h>

using namespace std;

void LocalSearch(int m, int n, vector<double> p,  list<int> locais_promissores, vector<double> &pen, vector<list<int>> &s, vector<vector<double>> c, vector<double> f, vector<double> d, vector<int> &atende, bool coeficientes_inteiros, double &otimo, double &UB, list<int> &abertas, list<int> &fechadas, double &fo, double soma_demandas, bool imprime_detalhes, bool &factivel, bool problema_esparso_e_gap, int Qmin, int R,  bool &resolveu_sp_na_otimalidade, int tempo_limite_heuristica, double &LB_modelo_otimo)
{
//cout << "******************************* APLICA HEURÍSTICA *******************************" << endl; //alterado

clock_t inicio_CPU_heuristica = clock();
int A = locais_promissores.size() ;
int B = locais_promissores.size()  ;
vector<int> indices_facilidades_abertas;
vector<int> vetor_locais_promissores;
for (std::list<int>::iterator k = locais_promissores.begin(); k != locais_promissores.end(); k++ ){
if (!s[*k].empty()  )
{
    indices_facilidades_abertas.push_back(*k);
    
}
vetor_locais_promissores.push_back(*k);
//indices_facilidades_abertas.push_back(*k);
}
list<int> clientes_que_foram_centro;
list<int> clientes_proibidos;
int Q = min(Qmin, A);
list<vector<int>> lista_tabu_clientes;
list<vector<int>> lista_tabu_facilidades; 
bool first_improvement = true;
bool best_improvement = false;
//int iter_sem_melhora = 0;

int MaiorQ = Qmin;
int IterSemMelhora = 0;
bool pare = false;
if (first_improvement == true)
{
while ( (IterSemMelhora <= R) && (Q <= A) && (tempo_decorrido(inicio_CPU_heuristica) <= tempo_limite_heuristica) )
{
bool pare = false;
int iter_sem_melhora = 0;
int iter2 = 0;
int ultima_melhora = 0;
vector<int> melhor_atende = atende;
double melhor_melhora = 0;
if (Q == A)
{
vector<int> vetor_facilidades_proximas1 = vetor_locais_promissores;  //vetor do cluster
vector<int> vetor_conjunto_clientes; 
    for (int j1 = 0; j1 < vetor_locais_promissores.size(); j1++)
    {
        for (std::list<int>::iterator k = s[vetor_facilidades_proximas1[j1]].begin(); k != s[vetor_facilidades_proximas1[j1]].end(); k++ )
        {
        //vetor_conjunto_clientes.push_back(*k);
        vetor_conjunto_clientes.push_back(*k);
    
        }
    }
bool factivel2 = true;
    for (int j = 0; j < vetor_facilidades_proximas1.size(); j++)
    {
       
       if (pen[vetor_facilidades_proximas1[j]] > 0 )
       {
           factivel2 = false;
           break;
       }
    }
 IloEnv env;
 IloModel mod(env);
 IloCplex cplex(mod);
 IloNumVarArray y(env, vetor_facilidades_proximas1.size(), 0, 1, ILOBOOL);
 IloArray<IloNumVarArray> x1(env, vetor_conjunto_clientes.size());
 //IloArray<IloNumVarArray> x2(env, vetor_conjunto_clientes.size());
for(int i = 0; i < vetor_conjunto_clientes.size(); i++) {
    x1[i] = IloNumVarArray(env, vetor_facilidades_proximas1.size(), 0, 1, ILOBOOL);
 }
 IloExpr expfo(env);
 for (int j = 0; j <  vetor_facilidades_proximas1.size(); j++){
	expfo += f[vetor_facilidades_proximas1[j]] * y[j];
 }
 for (int j = 0; j <  vetor_facilidades_proximas1.size(); j++){
		for(int i = 0; i < vetor_conjunto_clientes.size(); i++){
		expfo += c[vetor_conjunto_clientes[i]][vetor_facilidades_proximas1[j]]* d[vetor_conjunto_clientes[i]] * x1[i][j];
	}
 }
 IloNumVarArray h1(env, vetor_facilidades_proximas1.size(), 0, IloInfinity, ILOFLOAT);
if (factivel2 == false )
{
 for (int j = 0; j <  vetor_facilidades_proximas1.size(); j++){
	expfo += 100000*h1[j];
 }
}
            if (factivel2 == false)
            {
			for (int j = 0; j < vetor_facilidades_proximas1.size(); j++){
			IloExpr r3(env);
			for (int i = 0; i < vetor_conjunto_clientes.size(); i++){
				r3 += d[vetor_conjunto_clientes[i]]*x1[i][j];
			}
            r3 = r3 - h1[j];
			mod.add(r3  <= p[vetor_facilidades_proximas1[j]]*y[j]);  //tanto que cabe em cada uma
			r3.end();
		    }
            }
            else
            {
               for (int j = 0; j < vetor_facilidades_proximas1.size(); j++){
			IloExpr r3(env);
			for (int i = 0; i < vetor_conjunto_clientes.size(); i++){
				r3 += d[vetor_conjunto_clientes[i]]*x1[i][j];
			}
            //r3 = r3 - h1[j];
			mod.add(r3  <= p[vetor_facilidades_proximas1[j]]*y[j]);  //tanto que cabe em cada uma
			r3.end();
		    } 
            }

 IloAdd(mod, IloMinimize(env, expfo));
 expfo.end();

 //restricao de indivisibilidade*****************************************************************
 for (int i = 0; i < vetor_conjunto_clientes.size(); i++){
			IloExpr r2(env);
			for (int j = 0; j < vetor_facilidades_proximas1.size(); j++){
				r2 += x1[i][j];
			}
         	mod.add(r2 == 1);
			r2.end();
    }

if (problema_esparso_e_gap == true)
{
    IloExpr r8(env);
    for (int j = 0; j < vetor_facilidades_proximas1.size(); j++)
    {
        r8 += y[j];
    }
    mod.add(r8 == vetor_facilidades_proximas1.size());
    r8.end();
}
vector<vector<int>> start1(vetor_conjunto_clientes.size(), vector<int>(vetor_facilidades_proximas1.size()));
for (int i = 0; i < vetor_conjunto_clientes.size(); i++)
{
      for (int j = 0; j < vetor_facilidades_proximas1.size(); j++)
      {
          start1[i][j] = 0;
          if (atende[vetor_conjunto_clientes[i]] == vetor_facilidades_proximas1[j])
          {
              start1[i][j] = 1;

          }
      }
}
IloNumVarArray startVar(env);
IloNumArray startVal(env);
for (int i = 0; i < vetor_conjunto_clientes.size(); i++)
         for (int j = 0; j < vetor_facilidades_proximas1.size(); j++) {
             startVar.add(x1[i][j]);
             startVal.add(start1[i][j]);
         }
     cplex.addMIPStart(startVar, startVal);
     startVal.end();
     startVar.end();
cplex.setWarning(env.getNullStream()); // Eliminar warnings
cplex.setOut(env.getNullStream()); /// Eliminar os logs do solver
if (coeficientes_inteiros == true)
{
    cplex.setParam(IloCplex::EpGap, 0);
cplex.setParam(IloCplex::Param::MIP::Tolerances::AbsMIPGap, 0.99);
//cout << "coeficientes inteiros" << endl;
}

cplex.setParam(IloCplex::Param::Parallel,1);
cplex.setParam(IloCplex::Threads,1);
cplex.solve();



    for (int j1 = 0; j1 < vetor_facilidades_proximas1.size(); j1++)
    {
       for (int i = 0; i < vetor_conjunto_clientes.size(); i++)
     {
       if ( cplex.getValue(x1[i][j1]) > 0.9)
     {

       melhor_atende[vetor_conjunto_clientes[i]] = vetor_facilidades_proximas1[j1];
    }
    }
    }


atende = melhor_atende;
factivel = true;
abertas.clear();
fechadas.clear();
double fo1 = 0;
double fo2 = 0;
for (int j = 0; j < m; j++)
{
   s[j].clear();
   pen[j] = -p[j];
}
for (int i = 0; i < n; i++)
{
    s[atende[i]].push_back(i);
    pen[atende[i]] = pen[atende[i]] + d[i];
    fo1 = fo1 + d[i]*c[i][atende[i]]; 
}
for (int j = 0; j < m; j++)
{
   if (!s[j].empty())
   {
       fo2 = fo2 + f[j];
       abertas.push_back(j);
       if (pen[j] > 0)
       {
           fo2 = fo2 + pow(10,5)*pen[j];
           factivel = false;
       }
   }
   else
   {
       fechadas.push_back(j);
   }
   
}

//indices_facilidades_abertas = cria_vetor_copia_de_lista(abertas);
fo = fo1 + fo2;
cout << setprecision(12) <<  "Upper Bound : " << setprecision(20) <<  fo << endl;
LB_modelo_otimo = cplex.getBestObjValue();
UB = cplex.getObjValue();
//IterSemMelhora = 0;
resolveu_sp_na_otimalidade = true;
return;



}


start:
if ( tempo_decorrido(inicio_CPU_heuristica) + 10 > tempo_limite_heuristica) 
{
    if (factivel == true) 
    {
    return;
    }
    else
    {

     cout << "aplica novo procedimento:" << endl;
      //vector<int> indices_facilidades_abertas;

    funcao_ordena_decrescente(indices_facilidades_abertas,pen);

    double fo_antes = 0;

    double soma_capacidades_indices_facilidades_abertas = 0;
    for (int i = 0; i < indices_facilidades_abertas.size(); i++)
    {
        soma_capacidades_indices_facilidades_abertas = soma_capacidades_indices_facilidades_abertas + p[indices_facilidades_abertas[i]];
    }
    double soma_capacidades_locais_promissores = 0;
    for (int i = 0; i < vetor_locais_promissores.size(); i++)
    {
        soma_capacidades_locais_promissores = soma_capacidades_locais_promissores + p[vetor_locais_promissores[i]];
    }
    
    cout << "soma capacidades indices facilidades abertas: " << soma_capacidades_indices_facilidades_abertas << endl;
    cout << "soma capacidades vetor_locais_promissores: " << soma_capacidades_indices_facilidades_abertas << endl;
    cout << "soma demandas: " << soma_demandas << endl;
    cout << "indices_facilidades_abertas" << endl;
    for (int i = 0; i < indices_facilidades_abertas.size(); i++)
    {
        cout << indices_facilidades_abertas[i] << " ";
    }
    cout << endl;
        cout << "penalidades" << endl;
    for (int i = 0; i < indices_facilidades_abertas.size(); i++)
    {
        cout << pen[indices_facilidades_abertas[i]] << " ";
    }
    cout << endl;
    


/*
    imprime_vetor(atende);
    cout << endl;
    cout << "Vetor locais promissores" << endl;
    imprime_lista(locais_promissores);
    cout << endl;
*/
  
  
vector<int> indices_facilidades_abertas_violadas;
for (int i = 0; i < indices_facilidades_abertas.size(); i++)
{
    if (pen[indices_facilidades_abertas[i]] > 0) indices_facilidades_abertas_violadas.push_back(indices_facilidades_abertas[i]);
}

    IloEnv env;
    for (int j5 = 0; j5 < indices_facilidades_abertas_violadas.size() ; j5++)
    {

    //if (factivel == true) break;

    //if (pen[indices_facilidades_abertas_violadas[j5]] <= 0) continue;
    cout << "Facilidade:" << indices_facilidades_abertas_violadas[j5] << endl;     
// constroi o cluster
vector<int> vetor_facilidades_proximas1;  //vetor do cluster
vector<int> vetor_clientes_atendidos; ////novo
vector<int> vetor_conjunto_clientes_cluster;


    for (std::list<int>::iterator k = s[indices_facilidades_abertas_violadas[j5]].begin() ; k !=  s[indices_facilidades_abertas_violadas[j5]].end() ; k++)
    {
    vetor_clientes_atendidos.push_back(*k);
    }

    int cliente_mais_proximo;

//ordena clientes atendidos de acordo com a proximidade da facilidade centro
        for (int i = 0; i < vetor_clientes_atendidos.size(); i++)
    {
    int cliente_mais_proximo = i;
    for (int j = i; j < vetor_clientes_atendidos.size(); j++)
    {
    if (c[vetor_clientes_atendidos[i]][indices_facilidades_abertas_violadas[j5]]< c[vetor_clientes_atendidos[cliente_mais_proximo]][indices_facilidades_abertas_violadas[j5]])
        {
        cliente_mais_proximo = j;  //indice do mais proximo ocupa a posicao j
        }
    }
    //se mais proximo não está na posicao i
    if (cliente_mais_proximo != i)
    {
    int temp = vetor_clientes_atendidos[i];
    vetor_clientes_atendidos[i] = vetor_clientes_atendidos[cliente_mais_proximo];
    vetor_clientes_atendidos[cliente_mais_proximo] = temp;
    }
    }






    cliente_mais_proximo = vetor_clientes_atendidos[0]; 
    
   //ordena locais promissores na ordem de proximidade com relacao ao cliente mais proximo
for (int i = 0; i < vetor_locais_promissores.size(); i++)
{
    int facilidade_mais_proxima = i;
    for (int j = i; j < vetor_locais_promissores.size(); j++)
    {
    if (c[cliente_mais_proximo][vetor_locais_promissores[j]]< c[cliente_mais_proximo][vetor_locais_promissores[facilidade_mais_proxima]])
        {
        facilidade_mais_proxima = j;  //indice do mais proximo ocupa a posicao j
        }
    }
    //se mais proximo não está na posicao i
    if (facilidade_mais_proxima != i)
    {
    int temp = vetor_locais_promissores[i];
    vetor_locais_promissores[i] = vetor_locais_promissores[facilidade_mais_proxima];
    vetor_locais_promissores[facilidade_mais_proxima] = temp;
    }
 }


list<int> lista_facilidades_proximas1;
//vetor_facilidades_proximas1[0] = indices_facilidades_abertas_violadas[j5];
lista_facilidades_proximas1.push_back(indices_facilidades_abertas_violadas[j5]);
fo_antes = 0;
fo_antes = fo_antes +  f[indices_facilidades_abertas_violadas[j5]] + pow(10,5)*max(pen[indices_facilidades_abertas_violadas[j5]], 0.0);

int conta2 = 1;
int aux = 0;
int Q = Qmin;
//seleciona Q facilidades abertas mais proximas (lista de facilidades proximas)
while (conta2 < Q)
{
         if ( (vetor_locais_promissores[aux] != indices_facilidades_abertas_violadas[j5]) && (!s[ vetor_locais_promissores[aux]].empty()  ))
        {
        // if (indices_facilidades_abertas_violadas_copia[aux] != atende[i5]) //
        //{    //
            lista_facilidades_proximas1.push_back(vetor_locais_promissores[aux]);

            //vetor_facilidades_proximas1[conta2] = indices_facilidades_abertas_violadas_copia[aux];
            if (!s[vetor_locais_promissores[aux]].empty())
            {
            fo_antes = fo_antes + f[vetor_locais_promissores[aux]] + pow(10,5)*max(pen[vetor_locais_promissores[aux]],0.0)  ;
            //cout << indices_facilidades_abertas_violadas_copia[aux] << " esta aberta então adiciona o custo" << endl;   
            }
            conta2 = conta2 +1;
        }
    aux = aux + 1;
}

for (int j = 0; j < vetor_locais_promissores.size(); j++)
{

 if ((pen[vetor_locais_promissores[j]] < 0) && (!s[ vetor_locais_promissores[j]].empty()  ) && ( esta_na_lista( lista_facilidades_proximas1 , vetor_locais_promissores[j] ) == false ) ) 
 {
    lista_facilidades_proximas1.push_back(vetor_locais_promissores[j]);
    fo_antes = fo_antes + f[vetor_locais_promissores[j]] ;

 }
 if ( s[vetor_locais_promissores[j]].empty()  )  
 {
    lista_facilidades_proximas1.push_back(vetor_locais_promissores[j]);

 }
}





//cout << "vetor facilidades proximas 1: " << endl;
//imprime_vetor(vetor_facilidades_proximas1);
//cout << endl;
lista_facilidades_proximas1.sort();
vetor_facilidades_proximas1 = cria_vetor_copia_de_lista(lista_facilidades_proximas1);
list<int> lista_clientes;
for (int j = 0; j < vetor_facilidades_proximas1.size(); j++)
{
            for (std::list<int>::iterator k = s[vetor_facilidades_proximas1[j]].begin(); k != s[vetor_facilidades_proximas1[j]].end(); k++ )
        {
        //vetor_conjunto_clientes.push_back(*k);
        lista_clientes.push_back(*k);
        fo_antes = fo_antes + d[*k]*c[*k][vetor_facilidades_proximas1[j]];
        }
}

    lista_clientes.sort();
    vetor_conjunto_clientes_cluster = cria_vetor_copia_de_lista(lista_clientes);


bool factivel2 = true;
    for (int j = 0; j < vetor_facilidades_proximas1.size(); j++)
    {
       
       if (pen[vetor_facilidades_proximas1[j]] > 0 )
       {
           factivel2 = false;
           break;
       }
    }

vector<int> atende_copia(n);
for (int i = 0; i < n; i++)
{
    atende_copia[i] = atende[i];
}




 IloModel mod(env);
 IloCplex cplex(mod);
 IloNumVarArray y(env, vetor_facilidades_proximas1.size(), 0, 1, ILOBOOL);
 IloArray<IloNumVarArray> x1(env, vetor_conjunto_clientes_cluster.size());
 //IloArray<IloNumVarArray> x2(env, vetor_conjunto_clientes_cluster.size());

 for(int i = 0; i < vetor_conjunto_clientes_cluster.size(); i++) {
    x1[i] = IloNumVarArray(env, vetor_facilidades_proximas1.size(), 0, 1, ILOBOOL);
 }

 IloExpr expfo(env);
 for (int j = 0; j <  vetor_facilidades_proximas1.size(); j++){
	expfo += f[vetor_facilidades_proximas1[j]] * y[j];
 }
 for (int j = 0; j <  vetor_facilidades_proximas1.size(); j++){
		for(int i = 0; i < vetor_conjunto_clientes_cluster.size(); i++){
		expfo += c[vetor_conjunto_clientes_cluster[i]][vetor_facilidades_proximas1[j]]* d[vetor_conjunto_clientes_cluster[i]] * x1[i][j];
	}
 }

 IloNumVarArray h1(env, vetor_facilidades_proximas1.size(), 0, IloInfinity, ILOFLOAT);

if (factivel2 == false )
{
 for (int j = 0; j <  vetor_facilidades_proximas1.size(); j++){
	expfo += 100000*h1[j];
 }
}
            if (factivel2 == false)
            {
			for (int j = 0; j < vetor_facilidades_proximas1.size(); j++){
			IloExpr r3(env);
			for (int i = 0; i < vetor_conjunto_clientes_cluster.size(); i++){
				r3 += d[vetor_conjunto_clientes_cluster[i]]*x1[i][j];
			}
            r3 = r3 - h1[j];
			mod.add(r3  <= p[vetor_facilidades_proximas1[j]]*y[j]);  //tanto que cabe em cada uma
			r3.end();
		    }
            }
            else
            {
               for (int j = 0; j < vetor_facilidades_proximas1.size(); j++){
			IloExpr r3(env);
			for (int i = 0; i < vetor_conjunto_clientes_cluster.size(); i++){
				r3 += d[vetor_conjunto_clientes_cluster[i]]*x1[i][j];
			}
            //r3 = r3 - h1[j];
			mod.add(r3  <= p[vetor_facilidades_proximas1[j]]*y[j]);  //tanto que cabe em cada uma
			r3.end();
		    } 
            }

 IloAdd(mod, IloMinimize(env, expfo));
 expfo.end();

 //restricao de indivisibilidade*****************************************************************
 for (int i = 0; i < vetor_conjunto_clientes_cluster.size(); i++){
			IloExpr r2(env);
			for (int j = 0; j < vetor_facilidades_proximas1.size(); j++){
				r2 += x1[i][j];
			}
         	mod.add(r2 == 1);
			r2.end();
    }


/*
if (problema_esparso_e_gap == true)
{
    IloExpr r8(env);
    for (int j = 0; j < vetor_facilidades_proximas1.size(); j++)
    {
        r8 += y[j];
    }
    mod.add(r8 == vetor_facilidades_proximas1.size());
    r8.end();
}
*/



vector<vector<int>> start1(vetor_conjunto_clientes_cluster.size(), vector<int>(vetor_facilidades_proximas1.size()));
for (int i = 0; i < vetor_conjunto_clientes_cluster.size(); i++)
{
      for (int j = 0; j < vetor_facilidades_proximas1.size(); j++)
      {
          start1[i][j] = 0;
          if (atende[vetor_conjunto_clientes_cluster[i]] == vetor_facilidades_proximas1[j])
          {
              start1[i][j] = 1;

          }
      }
}

IloNumVarArray startVar(env);
IloNumArray startVal(env);
for (int i = 0; i < vetor_conjunto_clientes_cluster.size(); i++)
         for (int j = 0; j < vetor_facilidades_proximas1.size(); j++) {
             startVar.add(x1[i][j]);
             startVal.add(start1[i][j]);
         }
     cplex.addMIPStart(startVar, startVal);
     startVal.end();
     startVar.end();

//cplex.setWarning(env.getNullStream()); // Eliminar warnings
//cplex.setOut(env.getNullStream()); /// Eliminar os logs do solver
if (coeficientes_inteiros == true)
{
    cplex.setParam(IloCplex::EpGap, 0);
cplex.setParam(IloCplex::Param::MIP::Tolerances::AbsMIPGap, 0.99);
//cout << "coeficientes inteiros" << endl;
}



cout << "vetor facilidades proximas 1: " << vetor_facilidades_proximas1.size() << endl;
imprime_vetor(vetor_facilidades_proximas1);
cout << endl;
cout << "vetor conjunto clientes cluster: " << vetor_conjunto_clientes_cluster.size() <<  endl;
imprime_vetor(vetor_conjunto_clientes_cluster);
cout << endl;


cplex.setParam(IloCplex::Param::Parallel,1);
cplex.setParam(IloCplex::Threads,1);
//cout << "resolve subproblema.... Q = " << Q << " " << vetor_conjunto_clientes.size() << " x " << vetor_facilidades_proximas1.size() << " || "  << endl; //alterado

cplex.solve();


if(   (cplex.getStatus() == IloAlgorithm::Optimal) || (cplex.getStatus() == IloAlgorithm::Feasible) )
{
double fo_depois = cplex.getObjValue();

//cout << "fo_antes: " << fo_antes << endl;
//cout << "fo_depois: " << fo_depois << endl;

//cout << "fo_depois - fo_antes: " << fo_antes - fo_depois << endl;



if (fo_antes - fo_depois > 0.000001)  //houve melhora então atualiza solucao
{
 
 
     //cout << "atualiza solucao" << endl;
     //iter_sem_melhora = 0;
     //ultima_melhora = iter2;
     //clientes_proibidos.clear();
     //zera o clusteres repetidos

double melhor_melhora = fo_antes - fo_depois;
vector<int> melhor_atende = atende;
    for (int j1 = 0; j1 < vetor_facilidades_proximas1.size(); j1++)
    {
       for (int i = 0; i < vetor_conjunto_clientes_cluster.size(); i++)
     {
       if ( cplex.getValue(x1[i][j1]) > 0.9)
     {

       melhor_atende[vetor_conjunto_clientes_cluster[i]] = vetor_facilidades_proximas1[j1];
    }
    }
    }


atende = melhor_atende;
factivel = true;
abertas.clear();
fechadas.clear();
double fo1 = 0;
double fo2 = 0;
for (int j = 0; j < m; j++)
{
   s[j].clear();
   pen[j] = -p[j];
}
for (int i = 0; i < n; i++)
{
    s[atende[i]].push_back(i);
    pen[atende[i]] = pen[atende[i]] + d[i];
    fo1 = fo1 + d[i]*c[i][atende[i]]; 
}
for (int j = 0; j < m; j++)
{
   if (!s[j].empty())
   {
       fo2 = fo2 + f[j];
       abertas.push_back(j);
       if (pen[j] > 0)
       {
           fo2 = fo2 + pow(10,5)*pen[j];
           factivel = false;
       }
   }
   else
   {
       fechadas.push_back(j);
   }
   
}

//indices_facilidades_abertas = cria_vetor_copia_de_lista(abertas);
fo = fo1 + fo2;
cout << setprecision(12) <<  "Upper Bound : " << setprecision(20) <<  fo << endl;







}
}

mod.end();
cplex.end();



}

    return;


    }
}



cout << "Q: " << Q << endl;
bool houve_alguma_melhora = false; 
//cout << "recomeca" << endl;
vector<double> modulo_pen(m);

funcao_ordena_decrescente(indices_facilidades_abertas,pen);
vector<int> vetor_facilidades_proximas1;  //vetor do cluster
vector<int> vetor_conjunto_clientes; 
vector<int> vetor_clientes_atendidos; ////novo
IloEnv env;
for (int j5 = 0; j5 < indices_facilidades_abertas.size() ; j5++)
{

vetor_facilidades_proximas1.clear();
vetor_conjunto_clientes.clear();
vetor_clientes_atendidos.clear();
double fo_antes = 0;
bool prossegue = false;
bool continua = false;

int numero_tentativas = 0;
//cout << "faclidade: " << indices_facilidades_abertas[j5] << endl;
int cliente_mais_proximo;  //cliente qual será centralizado o cluster
bool conseguiu_montar_cluster = false; //novo 
if( !s[indices_facilidades_abertas[j5]].empty() )
    {

    for (std::list<int>::iterator k = s[indices_facilidades_abertas[j5]].begin() ; k !=  s[indices_facilidades_abertas[j5]].end() ; k++)
    {
    vetor_clientes_atendidos.push_back(*k);
    }      


    for (int i = 0; i < vetor_clientes_atendidos.size(); i++)
    {
    int cliente_mais_proximo = i;
    for (int j = i; j < vetor_clientes_atendidos.size(); j++)
    {
    if (c[vetor_clientes_atendidos[i]][indices_facilidades_abertas[j5]]< c[vetor_clientes_atendidos[cliente_mais_proximo]][indices_facilidades_abertas[j5]])
        {
        cliente_mais_proximo = j;  //indice do mais proximo ocupa a posicao j
        }
    }
    //se mais proximo não está na posicao i
    if (cliente_mais_proximo != i)
    {
    int temp = vetor_clientes_atendidos[i];
    vetor_clientes_atendidos[i] = vetor_clientes_atendidos[cliente_mais_proximo];
    vetor_clientes_atendidos[cliente_mais_proximo] = temp;
    }
    }

    
    }
  
int k = 0;    
//double fo_antes_auxiliar = fo_antes;
bool continuax = false;
while ( (conseguiu_montar_cluster == false) &&  (k < vetor_clientes_atendidos.size() )  )
{
    fo_antes = 0;
    vetor_facilidades_proximas1.clear();
    vetor_conjunto_clientes.clear();
    //cout << "valor de k: " << k << endl;
    //cout << "cliente mais proximo: " << vetor_clientes_atendidos[k] << endl; 
    conseguiu_montar_cluster = true;
    cliente_mais_proximo = vetor_clientes_atendidos[k]; 
    clientes_que_foram_centro.push_back(cliente_mais_proximo);
    clientes_proibidos.push_back(cliente_mais_proximo);
    
for (int i = 0; i < vetor_locais_promissores.size(); i++)
{
    int facilidade_mais_proxima = i;
    for (int j = i; j < vetor_locais_promissores.size(); j++)
    {
    if (c[cliente_mais_proximo][vetor_locais_promissores[j]]< c[cliente_mais_proximo][vetor_locais_promissores[facilidade_mais_proxima]])
        {
        facilidade_mais_proxima = j;  //indice do mais proximo ocupa a posicao j
        }
    }
    //se mais proximo não está na posicao i
    if (facilidade_mais_proxima != i)
    {
    int temp = vetor_locais_promissores[i];
    vetor_locais_promissores[i] = vetor_locais_promissores[facilidade_mais_proxima];
    vetor_locais_promissores[facilidade_mais_proxima] = temp;
    }
 }
list<int> lista_facilidades_proximas1;
//vetor_facilidades_proximas1[0] = indices_facilidades_abertas[j5];
lista_facilidades_proximas1.push_back(indices_facilidades_abertas[j5]);
fo_antes = fo_antes +  f[indices_facilidades_abertas[j5]] + pow(10,5)*max(pen[indices_facilidades_abertas[j5]], 0.0);

int conta2 = 1;
int aux = 0;

while (conta2 < Q)
{
         if ( (vetor_locais_promissores[aux] != indices_facilidades_abertas[j5]) && (!s[ vetor_locais_promissores[aux]].empty()  ))
        {
        // if (indices_facilidades_abertas_copia[aux] != atende[i5]) //
        //{    //
            lista_facilidades_proximas1.push_back(vetor_locais_promissores[aux]);

            //vetor_facilidades_proximas1[conta2] = indices_facilidades_abertas_copia[aux];
            if (!s[vetor_locais_promissores[aux]].empty())
            {
            fo_antes = fo_antes + f[vetor_locais_promissores[aux]] + pow(10,5)*max(pen[vetor_locais_promissores[aux]],0.0)  ;
            //cout << indices_facilidades_abertas_copia[aux] << " esta aberta então adiciona o custo" << endl;   
            }
            conta2 = conta2 +1;
        }
    aux = aux + 1;
}

lista_facilidades_proximas1.sort();
vetor_facilidades_proximas1 = cria_vetor_copia_de_lista(lista_facilidades_proximas1);
list<int> lista_clientes;
    for (int j1 = 0; j1 < Q; j1++)
    {
        for (std::list<int>::iterator k = s[vetor_facilidades_proximas1[j1]].begin(); k != s[vetor_facilidades_proximas1[j1]].end(); k++ )
        {
        //vetor_conjunto_clientes.push_back(*k);
        lista_clientes.push_back(*k);
        fo_antes = fo_antes + d[*k]*c[*k][vetor_facilidades_proximas1[j1]];
        }
    }

    lista_clientes.sort();
    vetor_conjunto_clientes = cria_vetor_copia_de_lista(lista_clientes);
bool facilidades_esta_na_lista_tabu = false;
bool clientes_esta_na_lista_tabu = false;
//verifica se vetor_conjunto_clientes_esta_na_lista
 auto it1 = std::find_if(lista_tabu_facilidades.begin(), lista_tabu_facilidades.end(),
                           [&](const std::vector<int>& v) { return v == vetor_facilidades_proximas1; });
    if (it1 != lista_tabu_facilidades.end()) {
        facilidades_esta_na_lista_tabu = true;
    } 
 auto it2 = std::find_if(lista_tabu_clientes.begin(), lista_tabu_clientes.end(),
                           [&](const std::vector<int>& v) { return v == vetor_conjunto_clientes; });
    if (it2 != lista_tabu_clientes.end()) {
        clientes_esta_na_lista_tabu = true;
    } 
if ((facilidades_esta_na_lista_tabu == true) && (clientes_esta_na_lista_tabu == true))
{
     //cout << "subproblema_repetido" << endl;
     iter_sem_melhora = iter_sem_melhora + 1;
     conseguiu_montar_cluster = false;
     k = 1000; //se k = 1000 não tenta novamente
     continuax = true;
}
  
k = k+1;
if (k > vetor_clientes_atendidos.size() )continuax = true;

} // fim while conseguiu montar cluster

if (continuax == true) continue;
bool factivel2 = true;
    for (int j = 0; j < vetor_facilidades_proximas1.size(); j++)
    {
       
       if (pen[vetor_facilidades_proximas1[j]] > 0 )
       {
           factivel2 = false;
           break;
       }
    }

vector<int> atende_copia(n);
for (int i = 0; i < n; i++) atende_copia[i] = atende[i];

 IloModel mod(env);
 IloCplex cplex(mod);
 IloNumVarArray y(env, vetor_facilidades_proximas1.size(), 0, 1, ILOBOOL);
 IloArray<IloNumVarArray> x1(env, vetor_conjunto_clientes.size());
 //IloArray<IloNumVarArray> x2(env, vetor_conjunto_clientes.size());

 for(int i = 0; i < vetor_conjunto_clientes.size(); i++) {
    x1[i] = IloNumVarArray(env, vetor_facilidades_proximas1.size(), 0, 1, ILOBOOL);
 }

 IloExpr expfo(env);
 for (int j = 0; j <  vetor_facilidades_proximas1.size(); j++){
	expfo += f[vetor_facilidades_proximas1[j]] * y[j];
 }
 for (int j = 0; j <  vetor_facilidades_proximas1.size(); j++){
		for(int i = 0; i < vetor_conjunto_clientes.size(); i++){
		expfo += c[vetor_conjunto_clientes[i]][vetor_facilidades_proximas1[j]]* d[vetor_conjunto_clientes[i]] * x1[i][j];
	}
 }

 IloNumVarArray h1(env, vetor_facilidades_proximas1.size(), 0, IloInfinity, ILOFLOAT);

if (factivel2 == false )
{
 for (int j = 0; j <  vetor_facilidades_proximas1.size(); j++){
	expfo += 100000*h1[j];
 }
}
            if (factivel2 == false)
            {
			for (int j = 0; j < vetor_facilidades_proximas1.size(); j++){
			IloExpr r3(env);
			for (int i = 0; i < vetor_conjunto_clientes.size(); i++){
				r3 += d[vetor_conjunto_clientes[i]]*x1[i][j];
			}
            r3 = r3 - h1[j];
			mod.add(r3  <= p[vetor_facilidades_proximas1[j]]*y[j]);  //tanto que cabe em cada uma
			r3.end();
		    }
            }
            else
            {
               for (int j = 0; j < vetor_facilidades_proximas1.size(); j++){
			IloExpr r3(env);
			for (int i = 0; i < vetor_conjunto_clientes.size(); i++){
				r3 += d[vetor_conjunto_clientes[i]]*x1[i][j];
			}
            //r3 = r3 - h1[j];
			mod.add(r3  <= p[vetor_facilidades_proximas1[j]]*y[j]);  //tanto que cabe em cada uma
			r3.end();
		    } 
            }

 IloAdd(mod, IloMinimize(env, expfo));
 expfo.end();

 //restricao de indivisibilidade*****************************************************************
 for (int i = 0; i < vetor_conjunto_clientes.size(); i++){
			IloExpr r2(env);
			for (int j = 0; j < vetor_facilidades_proximas1.size(); j++){
				r2 += x1[i][j];
			}
         	mod.add(r2 == 1);
			r2.end();
    }

if (problema_esparso_e_gap == true)
{
    IloExpr r8(env);
    for (int j = 0; j < vetor_facilidades_proximas1.size(); j++)
    {
        r8 += y[j];
    }
    mod.add(r8 == vetor_facilidades_proximas1.size());
    r8.end();
}


vector<vector<int>> start1(vetor_conjunto_clientes.size(), vector<int>(vetor_facilidades_proximas1.size()));
for (int i = 0; i < vetor_conjunto_clientes.size(); i++)
{
      for (int j = 0; j < vetor_facilidades_proximas1.size(); j++)
      {
          start1[i][j] = 0;
          if (atende[vetor_conjunto_clientes[i]] == vetor_facilidades_proximas1[j])
          {
              start1[i][j] = 1;

          }
      }
}

IloNumVarArray startVar(env);
IloNumArray startVal(env);
for (int i = 0; i < vetor_conjunto_clientes.size(); i++)
         for (int j = 0; j < vetor_facilidades_proximas1.size(); j++) {
             startVar.add(x1[i][j]);
             startVal.add(start1[i][j]);
         }
     cplex.addMIPStart(startVar, startVal);
     startVal.end();
     startVar.end();

cplex.setWarning(env.getNullStream()); // Eliminar warnings
cplex.setOut(env.getNullStream()); /// Eliminar os logs do solver
if (coeficientes_inteiros == true)
{
cplex.setParam(IloCplex::EpGap, 0);
cplex.setParam(IloCplex::Param::MIP::Tolerances::AbsMIPGap, 0.99);

//cout << "coeficientes inteiros" << endl;
}

//cplex.setParam(IloCplex::EpGap, 0.00001);


if ( tempo_decorrido(inicio_CPU_heuristica) + 10 > tempo_limite_heuristica)
{
    if (factivel == true) 
    {
        return;
    }
    else
    {



     cout << "Heuristic second phase: " << endl;
    funcao_ordena_decrescente(indices_facilidades_abertas,pen);
    double fo_antes = 0;

    double soma_capacidades_indices_facilidades_abertas = 0;
    for (int i = 0; i < indices_facilidades_abertas.size(); i++)
    {
        soma_capacidades_indices_facilidades_abertas = soma_capacidades_indices_facilidades_abertas + p[indices_facilidades_abertas[i]];
    }
    
vector<int> indices_facilidades_abertas_violadas;
for (int i = 0; i < indices_facilidades_abertas.size(); i++)
{
    if (pen[indices_facilidades_abertas[i]] > 0) indices_facilidades_abertas_violadas.push_back(indices_facilidades_abertas[i]);
}

    IloEnv env;
    for (int j5 = 0; j5 < indices_facilidades_abertas_violadas.size() ; j5++)
    {
vector<int> vetor_facilidades_proximas1;  //vetor do cluster
vector<int> vetor_clientes_atendidos; ////novo
vector<int> vetor_conjunto_clientes_cluster;
    for (std::list<int>::iterator k = s[indices_facilidades_abertas_violadas[j5]].begin() ; k !=  s[indices_facilidades_abertas_violadas[j5]].end() ; k++)
    {
    vetor_clientes_atendidos.push_back(*k);
    }
    int cliente_mais_proximo;
        for (int i = 0; i < vetor_clientes_atendidos.size(); i++)
    {
    int cliente_mais_proximo = i;
    for (int j = i; j < vetor_clientes_atendidos.size(); j++)
    {
    if (c[vetor_clientes_atendidos[i]][indices_facilidades_abertas_violadas[j5]]< c[vetor_clientes_atendidos[cliente_mais_proximo]][indices_facilidades_abertas_violadas[j5]])
        {
        cliente_mais_proximo = j;  //indice do mais proximo ocupa a posicao j
        }
    }
    if (cliente_mais_proximo != i)
    {
    int temp = vetor_clientes_atendidos[i];
    vetor_clientes_atendidos[i] = vetor_clientes_atendidos[cliente_mais_proximo];
    vetor_clientes_atendidos[cliente_mais_proximo] = temp;
    }
    }
    cliente_mais_proximo = vetor_clientes_atendidos[0]; 
for (int i = 0; i < vetor_locais_promissores.size(); i++)
{
    int facilidade_mais_proxima = i;
    for (int j = i; j < vetor_locais_promissores.size(); j++)
    {
    if (c[cliente_mais_proximo][vetor_locais_promissores[j]]< c[cliente_mais_proximo][vetor_locais_promissores[facilidade_mais_proxima]])
        {
        facilidade_mais_proxima = j;  //indice do mais proximo ocupa a posicao j
        }
    }
    if (facilidade_mais_proxima != i)
    {
    int temp = vetor_locais_promissores[i];
    vetor_locais_promissores[i] = vetor_locais_promissores[facilidade_mais_proxima];
    vetor_locais_promissores[facilidade_mais_proxima] = temp;
    }
 }
list<int> lista_facilidades_proximas1;
lista_facilidades_proximas1.push_back(indices_facilidades_abertas_violadas[j5]);
fo_antes = 0;
fo_antes = fo_antes +  f[indices_facilidades_abertas_violadas[j5]] + pow(10,5)*max(pen[indices_facilidades_abertas_violadas[j5]], 0.0);
int conta2 = 1;
int aux = 0;
int Q = Qmin;
while (conta2 < Q)
{
         if ( (vetor_locais_promissores[aux] != indices_facilidades_abertas_violadas[j5]) && (!s[ vetor_locais_promissores[aux]].empty()  ))
        {
            lista_facilidades_proximas1.push_back(vetor_locais_promissores[aux]);
            if (!s[vetor_locais_promissores[aux]].empty())
            {
            fo_antes = fo_antes + f[vetor_locais_promissores[aux]] + pow(10,5)*max(pen[vetor_locais_promissores[aux]],0.0)  ;
            }
            conta2 = conta2 +1;
        }
    aux = aux + 1;
}
for (int j = 0; j < vetor_locais_promissores.size(); j++)
{

 if ((pen[vetor_locais_promissores[j]] < 0) && (!s[ vetor_locais_promissores[j]].empty()  ) && ( esta_na_lista( lista_facilidades_proximas1 , vetor_locais_promissores[j] ) == false ) ) 
 {
    lista_facilidades_proximas1.push_back(vetor_locais_promissores[j]);
   fo_antes = fo_antes + f[vetor_locais_promissores[j]] ;
 }  
}
lista_facilidades_proximas1.sort();
vetor_facilidades_proximas1 = cria_vetor_copia_de_lista(lista_facilidades_proximas1);
list<int> lista_clientes;
for (int j = 0; j < vetor_facilidades_proximas1.size(); j++)
{
            for (std::list<int>::iterator k = s[vetor_facilidades_proximas1[j]].begin(); k != s[vetor_facilidades_proximas1[j]].end(); k++ )
        {
        lista_clientes.push_back(*k);
        fo_antes = fo_antes + d[*k]*c[*k][vetor_facilidades_proximas1[j]];
        }
}
    lista_clientes.sort();
    vetor_conjunto_clientes_cluster = cria_vetor_copia_de_lista(lista_clientes);
bool factivel2 = true;
    for (int j = 0; j < vetor_facilidades_proximas1.size(); j++)
    {
       if (pen[vetor_facilidades_proximas1[j]] > 0 )
       {
           factivel2 = false;
           break;
       }
    }
vector<int> atende_copia(n);
for (int i = 0; i < n; i++)
{
    atende_copia[i] = atende[i];
}
 IloModel mod(env);
 IloCplex cplex(mod);
 IloNumVarArray y(env, vetor_facilidades_proximas1.size(), 0, 1, ILOBOOL);
 IloArray<IloNumVarArray> x1(env, vetor_conjunto_clientes_cluster.size());
 for(int i = 0; i < vetor_conjunto_clientes_cluster.size(); i++) {
    x1[i] = IloNumVarArray(env, vetor_facilidades_proximas1.size(), 0, 1, ILOBOOL);
 }
 IloExpr expfo(env);
 for (int j = 0; j <  vetor_facilidades_proximas1.size(); j++){
	expfo += f[vetor_facilidades_proximas1[j]] * y[j];
 }
 for (int j = 0; j <  vetor_facilidades_proximas1.size(); j++){
		for(int i = 0; i < vetor_conjunto_clientes_cluster.size(); i++){
		expfo += c[vetor_conjunto_clientes_cluster[i]][vetor_facilidades_proximas1[j]]* d[vetor_conjunto_clientes_cluster[i]] * x1[i][j];
	}
 }
 IloNumVarArray h1(env, vetor_facilidades_proximas1.size(), 0, IloInfinity, ILOFLOAT);
if (factivel2 == false )
{
 for (int j = 0; j <  vetor_facilidades_proximas1.size(); j++){
	expfo += 100000*h1[j];
 }
}
            if (factivel2 == false)
            {
			for (int j = 0; j < vetor_facilidades_proximas1.size(); j++){
			IloExpr r3(env);
			for (int i = 0; i < vetor_conjunto_clientes_cluster.size(); i++){
				r3 += d[vetor_conjunto_clientes_cluster[i]]*x1[i][j];
			}
            r3 = r3 - h1[j];
			mod.add(r3  <= p[vetor_facilidades_proximas1[j]]*y[j]);  //tanto que cabe em cada uma
			r3.end();
		    }
            }
            else
            {
               for (int j = 0; j < vetor_facilidades_proximas1.size(); j++){
			IloExpr r3(env);
			for (int i = 0; i < vetor_conjunto_clientes_cluster.size(); i++){
				r3 += d[vetor_conjunto_clientes_cluster[i]]*x1[i][j];
			}
            //r3 = r3 - h1[j];
			mod.add(r3  <= p[vetor_facilidades_proximas1[j]]*y[j]);  //tanto que cabe em cada uma
			r3.end();
		    } 
            }

 IloAdd(mod, IloMinimize(env, expfo));
 expfo.end();
 for (int i = 0; i < vetor_conjunto_clientes_cluster.size(); i++){
			IloExpr r2(env);
			for (int j = 0; j < vetor_facilidades_proximas1.size(); j++){
				r2 += x1[i][j];
			}
         	mod.add(r2 == 1);
			r2.end();
    }
if (problema_esparso_e_gap == true)
{
    IloExpr r8(env);
    for (int j = 0; j < vetor_facilidades_proximas1.size(); j++)
    {
        r8 += y[j];
    }
    mod.add(r8 == vetor_facilidades_proximas1.size());
    r8.end();
}
vector<vector<int>> start1(vetor_conjunto_clientes_cluster.size(), vector<int>(vetor_facilidades_proximas1.size()));
for (int i = 0; i < vetor_conjunto_clientes_cluster.size(); i++)
{
      for (int j = 0; j < vetor_facilidades_proximas1.size(); j++)
      {
          start1[i][j] = 0;
          if (atende[vetor_conjunto_clientes_cluster[i]] == vetor_facilidades_proximas1[j])
          {
              start1[i][j] = 1;

          }
      }
}
IloNumVarArray startVar(env);
IloNumArray startVal(env);
for (int i = 0; i < vetor_conjunto_clientes_cluster.size(); i++)
         for (int j = 0; j < vetor_facilidades_proximas1.size(); j++) {
             startVar.add(x1[i][j]);
             startVal.add(start1[i][j]);
         }
     cplex.addMIPStart(startVar, startVal);
     startVal.end();
     startVar.end();
if (coeficientes_inteiros == true)
{
    cplex.setParam(IloCplex::EpGap, 0);
cplex.setParam(IloCplex::Param::MIP::Tolerances::AbsMIPGap, 0.99);
}
cplex.setParam(IloCplex::Param::Parallel,1);
cplex.setParam(IloCplex::Threads,1);
cplex.solve();
if(   (cplex.getStatus() == IloAlgorithm::Optimal) || (cplex.getStatus() == IloAlgorithm::Feasible) )
{
double fo_depois = cplex.getObjValue();
if (fo_antes - fo_depois > 0.000001)  //houve melhora então atualiza solucao
{
 
double melhor_melhora = fo_antes - fo_depois;
vector<int> melhor_atende = atende;
    for (int j1 = 0; j1 < vetor_facilidades_proximas1.size(); j1++)
    {
       for (int i = 0; i < vetor_conjunto_clientes_cluster.size(); i++)
     {
       if ( cplex.getValue(x1[i][j1]) > 0.9)
     {
       melhor_atende[vetor_conjunto_clientes_cluster[i]] = vetor_facilidades_proximas1[j1];
    }
    }
    }
atende = melhor_atende;
factivel = true;
abertas.clear();
fechadas.clear();
double fo1 = 0;
double fo2 = 0;
for (int j = 0; j < m; j++)
{
   s[j].clear();
   pen[j] = -p[j];
}
for (int i = 0; i < n; i++)
{
    s[atende[i]].push_back(i);
    pen[atende[i]] = pen[atende[i]] + d[i];
    fo1 = fo1 + d[i]*c[i][atende[i]]; 
}
for (int j = 0; j < m; j++)
{
   if (!s[j].empty())
   {
       fo2 = fo2 + f[j];
       abertas.push_back(j);
       if (pen[j] > 0)
       {
           fo2 = fo2 + pow(10,5)*pen[j];
           factivel = false;
       }
   }
   else
   {
       fechadas.push_back(j);
   }
   
}
fo = fo1 + fo2;
cout << setprecision(12) <<  "Upper Bound : " << setprecision(20) <<  fo << endl;
}
}
mod.end();
cplex.end();
}
    return;
    }
} 
double tilim = tempo_limite_heuristica - tempo_decorrido(inicio_CPU_heuristica);
cplex.setParam(IloCplex::TiLim, tilim);
cplex.setParam(IloCplex::Param::Parallel,1);
cplex.setParam(IloCplex::Threads,1);
clock_t start_CPU = clock();
cplex.solve();
clock_t end_CPU = clock();
lista_tabu_clientes.push_back(vetor_conjunto_clientes);
lista_tabu_facilidades.push_back(vetor_facilidades_proximas1);
if(   (cplex.getStatus() == IloAlgorithm::Optimal) || (cplex.getStatus() == IloAlgorithm::Feasible) )
{
double fo_depois = cplex.getObjValue();
if (fo_antes - fo_depois > 0.000001)  //houve melhora então atualiza solucao
{
melhor_melhora = fo_antes - fo_depois;
melhor_atende = atende;
    for (int j1 = 0; j1 < vetor_facilidades_proximas1.size(); j1++)
    {
       for (int i = 0; i < vetor_conjunto_clientes.size(); i++)
     {
       if ( cplex.getValue(x1[i][j1]) > 0.9)
     {
       melhor_atende[vetor_conjunto_clientes[i]] = vetor_facilidades_proximas1[j1];
    }
    }
    }
atende = melhor_atende;
factivel = true;
abertas.clear();
fechadas.clear();
double fo1 = 0;
double fo2 = 0;
for (int j = 0; j < m; j++)
{
   s[j].clear();
   pen[j] = -p[j];
}
for (int i = 0; i < n; i++)
{
    s[atende[i]].push_back(i);
    pen[atende[i]] = pen[atende[i]] + d[i];
    fo1 = fo1 + d[i]*c[i][atende[i]]; 
}
for (int j = 0; j < m; j++)
{
   if (!s[j].empty())
   {
       fo2 = fo2 + f[j];
       abertas.push_back(j);
       if (pen[j] > 0)
       {
           fo2 = fo2 + pow(10,5)*pen[j];
           factivel = false;
       }
   }
   else
   {
       fechadas.push_back(j);
   }
   
}

indices_facilidades_abertas = cria_vetor_copia_de_lista(abertas);
fo = fo1 + fo2;
cout << setprecision(12) <<  "Upper Bound : " << setprecision(20) <<  fo << endl;
IterSemMelhora = 0;

if (Q > MaiorQ)
{
MaiorQ = Q;
}

Q = Qmin; 
houve_alguma_melhora = true;

//env.end();
mod.end();
cplex.end();

goto start;
}
}
mod.end();
cplex.end();
}
env.end();
if (houve_alguma_melhora == false)
{
    
    if (Q > MaiorQ)
    {
    IterSemMelhora = IterSemMelhora + 1;
    }
    Q = Q+1;
    
}
}
}
if (factivel == false)
{
   cout << "Heuristic second phase:" << endl;
    funcao_ordena_decrescente(indices_facilidades_abertas,pen);
    double fo_antes = 0;
    double soma_capacidades_indices_facilidades_abertas = 0;
    for (int i = 0; i < indices_facilidades_abertas.size(); i++)
    {
        soma_capacidades_indices_facilidades_abertas = soma_capacidades_indices_facilidades_abertas + p[indices_facilidades_abertas[i]];
    }
vector<int> indices_facilidades_abertas_violadas;
for (int i = 0; i < indices_facilidades_abertas.size(); i++)
{
    if (pen[indices_facilidades_abertas[i]] > 0) indices_facilidades_abertas_violadas.push_back(indices_facilidades_abertas[i]);
}
    IloEnv env;
    for (int j5 = 0; j5 < indices_facilidades_abertas_violadas.size() ; j5++)
    {
vector<int> vetor_facilidades_proximas1;  //vetor do cluster
vector<int> vetor_clientes_atendidos; ////novo
vector<int> vetor_conjunto_clientes_cluster;
    for (std::list<int>::iterator k = s[indices_facilidades_abertas_violadas[j5]].begin() ; k !=  s[indices_facilidades_abertas_violadas[j5]].end() ; k++)
    {
    vetor_clientes_atendidos.push_back(*k);
    }
    int cliente_mais_proximo;
        for (int i = 0; i < vetor_clientes_atendidos.size(); i++)
    {
    int cliente_mais_proximo = i;
    for (int j = i; j < vetor_clientes_atendidos.size(); j++)
    {
    if (c[vetor_clientes_atendidos[i]][indices_facilidades_abertas_violadas[j5]]< c[vetor_clientes_atendidos[cliente_mais_proximo]][indices_facilidades_abertas_violadas[j5]])
        {
        cliente_mais_proximo = j;  //indice do mais proximo ocupa a posicao j
        }
    }
    //se mais proximo não está na posicao i
    if (cliente_mais_proximo != i)
    {
    int temp = vetor_clientes_atendidos[i];
    vetor_clientes_atendidos[i] = vetor_clientes_atendidos[cliente_mais_proximo];
    vetor_clientes_atendidos[cliente_mais_proximo] = temp;
    }
    }
    cliente_mais_proximo = vetor_clientes_atendidos[0]; 
for (int i = 0; i < vetor_locais_promissores.size(); i++)
{
    int facilidade_mais_proxima = i;
    for (int j = i; j < vetor_locais_promissores.size(); j++)
    {
    if (c[cliente_mais_proximo][vetor_locais_promissores[j]]< c[cliente_mais_proximo][vetor_locais_promissores[facilidade_mais_proxima]])
        {
        facilidade_mais_proxima = j;  //indice do mais proximo ocupa a posicao j
        }
    }
    if (facilidade_mais_proxima != i)
    {
    int temp = vetor_locais_promissores[i];
    vetor_locais_promissores[i] = vetor_locais_promissores[facilidade_mais_proxima];
    vetor_locais_promissores[facilidade_mais_proxima] = temp;
    }
 }
list<int> lista_facilidades_proximas1;
lista_facilidades_proximas1.push_back(indices_facilidades_abertas_violadas[j5]);
fo_antes = 0;
fo_antes = fo_antes +  f[indices_facilidades_abertas_violadas[j5]] + pow(10,5)*max(pen[indices_facilidades_abertas_violadas[j5]], 0.0);
int conta2 = 1;
int aux = 0;
int Q = Qmin;
while (conta2 < Q)
{
         if ( (vetor_locais_promissores[aux] != indices_facilidades_abertas_violadas[j5]) && (!s[ vetor_locais_promissores[aux]].empty()  ))
        {
            lista_facilidades_proximas1.push_back(vetor_locais_promissores[aux]);
            if (!s[vetor_locais_promissores[aux]].empty())
            {
            fo_antes = fo_antes + f[vetor_locais_promissores[aux]] + pow(10,5)*max(pen[vetor_locais_promissores[aux]],0.0)  ;
            }
            conta2 = conta2 +1;
        }
    aux = aux + 1;
}
for (int j = 0; j < vetor_locais_promissores.size(); j++)
{
 if ((pen[vetor_locais_promissores[j]] < 0) && (!s[ vetor_locais_promissores[j]].empty()  ) && ( esta_na_lista( lista_facilidades_proximas1 , vetor_locais_promissores[j] ) == false ) ) 
 {
    lista_facilidades_proximas1.push_back(vetor_locais_promissores[j]);
    fo_antes = fo_antes + f[vetor_locais_promissores[j]] ;
 }  
}
lista_facilidades_proximas1.sort();
vetor_facilidades_proximas1 = cria_vetor_copia_de_lista(lista_facilidades_proximas1);
list<int> lista_clientes;
for (int j = 0; j < vetor_facilidades_proximas1.size(); j++)
{
            for (std::list<int>::iterator k = s[vetor_facilidades_proximas1[j]].begin(); k != s[vetor_facilidades_proximas1[j]].end(); k++ )
        {
        lista_clientes.push_back(*k);
        fo_antes = fo_antes + d[*k]*c[*k][vetor_facilidades_proximas1[j]];
        }
}
    lista_clientes.sort();
    vetor_conjunto_clientes_cluster = cria_vetor_copia_de_lista(lista_clientes);
bool factivel2 = true;
    for (int j = 0; j < vetor_facilidades_proximas1.size(); j++)
    {
       if (pen[vetor_facilidades_proximas1[j]] > 0 )
       {
           factivel2 = false;
           break;
       }
    }
vector<int> atende_copia(n);
for (int i = 0; i < n; i++)
{
    atende_copia[i] = atende[i];
}
 IloModel mod(env);
 IloCplex cplex(mod);
 IloNumVarArray y(env, vetor_facilidades_proximas1.size(), 0, 1, ILOBOOL);
 IloArray<IloNumVarArray> x1(env, vetor_conjunto_clientes_cluster.size());
 for(int i = 0; i < vetor_conjunto_clientes_cluster.size(); i++) {
    x1[i] = IloNumVarArray(env, vetor_facilidades_proximas1.size(), 0, 1, ILOBOOL);
 }
 IloExpr expfo(env);
 for (int j = 0; j <  vetor_facilidades_proximas1.size(); j++){
	expfo += f[vetor_facilidades_proximas1[j]] * y[j];
 }
 for (int j = 0; j <  vetor_facilidades_proximas1.size(); j++){
		for(int i = 0; i < vetor_conjunto_clientes_cluster.size(); i++){
		expfo += c[vetor_conjunto_clientes_cluster[i]][vetor_facilidades_proximas1[j]]* d[vetor_conjunto_clientes_cluster[i]] * x1[i][j];
	}
 }
 IloNumVarArray h1(env, vetor_facilidades_proximas1.size(), 0, IloInfinity, ILOFLOAT);
if (factivel2 == false )
{
 for (int j = 0; j <  vetor_facilidades_proximas1.size(); j++){
	expfo += 100000*h1[j];
 }
}
            if (factivel2 == false)
            {
			for (int j = 0; j < vetor_facilidades_proximas1.size(); j++){
			IloExpr r3(env);
			for (int i = 0; i < vetor_conjunto_clientes_cluster.size(); i++){
				r3 += d[vetor_conjunto_clientes_cluster[i]]*x1[i][j];
			}
            r3 = r3 - h1[j];
			mod.add(r3  <= p[vetor_facilidades_proximas1[j]]*y[j]);  //tanto que cabe em cada uma
			r3.end();
		    }
            }
            else
            {
               for (int j = 0; j < vetor_facilidades_proximas1.size(); j++){
			IloExpr r3(env);
			for (int i = 0; i < vetor_conjunto_clientes_cluster.size(); i++){
				r3 += d[vetor_conjunto_clientes_cluster[i]]*x1[i][j];
			}
            //r3 = r3 - h1[j];
			mod.add(r3  <= p[vetor_facilidades_proximas1[j]]*y[j]);  //tanto que cabe em cada uma
			r3.end();
		    } 
            }

 IloAdd(mod, IloMinimize(env, expfo));
 expfo.end();
 for (int i = 0; i < vetor_conjunto_clientes_cluster.size(); i++){
			IloExpr r2(env);
			for (int j = 0; j < vetor_facilidades_proximas1.size(); j++){
				r2 += x1[i][j];
			}
         	mod.add(r2 == 1);
			r2.end();
    }
if (problema_esparso_e_gap == true)
{
    IloExpr r8(env);
    for (int j = 0; j < vetor_facilidades_proximas1.size(); j++)
    {
        r8 += y[j];
    }
    mod.add(r8 == vetor_facilidades_proximas1.size());
    r8.end();
}
vector<vector<int>> start1(vetor_conjunto_clientes_cluster.size(), vector<int>(vetor_facilidades_proximas1.size()));
for (int i = 0; i < vetor_conjunto_clientes_cluster.size(); i++)
{
      for (int j = 0; j < vetor_facilidades_proximas1.size(); j++)
      {
          start1[i][j] = 0;
          if (atende[vetor_conjunto_clientes_cluster[i]] == vetor_facilidades_proximas1[j])
          {
              start1[i][j] = 1;

          }
      }
}
IloNumVarArray startVar(env);
IloNumArray startVal(env);
for (int i = 0; i < vetor_conjunto_clientes_cluster.size(); i++)
         for (int j = 0; j < vetor_facilidades_proximas1.size(); j++) {
             startVar.add(x1[i][j]);
             startVal.add(start1[i][j]);
         }
     cplex.addMIPStart(startVar, startVal);
     startVal.end();
     startVar.end();

if (coeficientes_inteiros == true)
{
    cplex.setParam(IloCplex::EpGap, 0);
cplex.setParam(IloCplex::Param::MIP::Tolerances::AbsMIPGap, 0.99);
}

cplex.setParam(IloCplex::Param::Parallel,1);
cplex.setParam(IloCplex::Threads,1);
cplex.solve();
if(   (cplex.getStatus() == IloAlgorithm::Optimal) || (cplex.getStatus() == IloAlgorithm::Feasible) )
{
double fo_depois = cplex.getObjValue();
if (fo_antes - fo_depois > 0.000001)  //houve melhora então atualiza solucao
{
double melhor_melhora = fo_antes - fo_depois;
vector<int> melhor_atende = atende;
    for (int j1 = 0; j1 < vetor_facilidades_proximas1.size(); j1++)
    {
       for (int i = 0; i < vetor_conjunto_clientes_cluster.size(); i++)
     {
       if ( cplex.getValue(x1[i][j1]) > 0.9)
     {
       melhor_atende[vetor_conjunto_clientes_cluster[i]] = vetor_facilidades_proximas1[j1];
    }
    }
    }
atende = melhor_atende;
factivel = true;
abertas.clear();
fechadas.clear();
double fo1 = 0;
double fo2 = 0;
for (int j = 0; j < m; j++)
{
   s[j].clear();
   pen[j] = -p[j];
}
for (int i = 0; i < n; i++)
{
    s[atende[i]].push_back(i);
    pen[atende[i]] = pen[atende[i]] + d[i];
    fo1 = fo1 + d[i]*c[i][atende[i]]; 
}
for (int j = 0; j < m; j++)
{
   if (!s[j].empty())
   {
       fo2 = fo2 + f[j];
       abertas.push_back(j);
       if (pen[j] > 0)
       {
           fo2 = fo2 + pow(10,5)*pen[j];
           factivel = false;
       }
   }
   else
   {
       fechadas.push_back(j);
   }
   
}
fo = fo1 + fo2;
cout << setprecision(12) <<  "Upper Bound : " << setprecision(20) <<  fo << endl;
}
}
mod.end();
cplex.end();
}
    return;
}
}

