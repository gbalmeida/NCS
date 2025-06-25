// arquivo funcoes.h

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
#include <ctime>

#include <ilcplex/ilocplex.h>

using namespace std;

void Read_instances(char *argv, int &m, int &n, vector<double> &p, vector<double> &f, vector<double> &d,  vector<vector<double>> &c,  vector<double> &pen, double &demanda_total, double &capacidade_total, bool &coeficientes_inteiros) ;


double tempo_decorrido(clock_t& inicio_CPU);
void SolveSparseProblem(bool problema_esparso_e_gap, bool solucao_heuristica_e_factivel, int quantidade_iteracoes , vector<int> vetor_locais_promissores, vector<double> f, vector<int> vetor_conjunto_de_clientes, vector<vector<double>> c, vector<double> p, vector<double> pen, vector<double> d, double &UB, bool coeficientes_inteiros, int tempo_limite_problema_esparso, int &iteracao_que_encontrou_otimo,  bool &solucao_heuristica, vector<int> atende, list<vector<double>> &coeficientes_cortes_x, list<vector<int>> &clientes_cortes_x, vector<double> &limitantes_cortes_x, vector<int> &variavel_y_associada_ao_corte, list<vector<double>> &coeficientes_cortes_y, list<vector<int>> &facilidades_cortes_y, vector<double> &limitantes_cortes_y, vector<int> &melhor_atende, int n, int m, int tamanho_minimo, double &lb_esparso, double &LB_modelo_otimo);
void monta_e_resolve_lp();
bool verifica_corte_e_valido_usando_solver(vector<double> Pesos , double B, vector<double> coef, vector<int> clientes_corte, double limitante);
void aplica_fenchel_cut_clientes_novo (vector<double> Pesos , vector<int> Itens, double B, vector<double> &coef, vector<int> &clientes_corte, double &limitante, vector<double> &x, bool &corte_violado, bool &impossivel_cortar);
void FCColumnGeneration (vector<double> Pesos , vector<int> Itens, double B, vector<double> &coef, vector<int> &clientes_corte, double &limitante, vector<double> &x, bool &corte_violado, bool &impossivel_cortar, int &tentativas_opcao0, int &total_tentativas, int bigM, int max_iter_fc);
void CIgreedyA(vector<int> vetor_N1, vector<int> vetor_N2, vector<double> aux_denso, double &soma_demanda_cobertura, double B, vector<double> d, vector<int> &vetor_cobertura , bool &corte_encontrado, vector<double> x, double y_til, vector<int> &best_CI, double &valor_best_CI    );
void CIgreedyB(vector<int> vetor_N1, vector<int> vetor_N2, vector<double> aux_denso, double &soma_demanda_cobertura, double B, vector<double> d, vector<int> &vetor_cobertura, bool &corte_encontrado, vector<double> x, double y_til,  vector<int> &best_CI, double &valor_best_CI);
void AddOneCostumer(vector<int> vetor_N1, vector<int> vetor_N2, vector<double> aux_denso, double &soma_demanda_cobertura, double B, vector<double> d, vector<int> &vetor_cobertura, bool &corte_encontrado, vector<double> x ,double y_til,  vector<int> &best_CI, double &valor_best_CI  );
void AddTwoCostumers(vector<int> vetor_N1, vector<int> vetor_N2, vector<double> aux_denso, double &soma_demanda_cobertura, double B, vector<double> d, vector<int> &vetor_cobertura, bool &corte_encontrado, vector<double> x, double y_til,  vector<int> &best_CI, double &valor_best_CI   );
void AddThreeCostumers(vector<int> vetor_N1, vector<int> vetor_N2, vector<double> aux_denso, double &soma_demanda_cobertura, double B, vector<double> d, vector<int> &vetor_cobertura, bool &corte_encontrado, vector<double> x, double y_til,  vector<int> &best_CI, double &valor_best_CI   );
void CIexact(double B, vector<int> vetor_N1, vector<int> vetor_N2, vector<int> vetor_D , double y_til, vector<double> x, vector<double> d, vector<int> &vetor_cobertura_auxiliar, double &soma_demandas_cobertura_auxiliar, bool &corte_encontrado, double soma_demandas_D, bool &cobertura_encontrada,  vector<int> &best_CI, double &valor_best_CI);
bool esta_na_lista_double(const std::list<double>& lista, int numero);   //verifica se um número está em uma lista
void LocalSearch(int m, int n, vector<double> p,  list<int> locais_promissores, vector<double> &pen, vector<list<int>> &s, vector<vector<double>> c, vector<double> f, vector<double> d, vector<int> &atende, bool coeficientes_inteiros, double &otimo, double &UB, list<int> &abertas, list<int> &fechadas, double &fo, double soma_demandas, bool imprime_detalhes, bool &factivel, bool problema_esparso_e_gap, int Qmin, int R, bool &resolveu_sp_na_otimalidade, int tempo_limite_heuristica, double &LB_modelo_otimo);
int knapsack(int Capacidade, vector<int>& Itens, vector<int>& Pesos, vector<int>& Valores);
int mochila(double capacidade, vector<double> coef, vector<int> clientes_corte , vector<double> Pesos);
void SolveKnapsack(double capacidade, vector<double> coef, vector<int> clientes_corte , vector<double> Pesos, double &otimo_mochila, vector<int> &itens_otimo);
void imprime_lista_double(const std::list<double>& minhaLista) ;

void LCIsuperaditivefunction(vector<double> Pesos , vector<int> Itens, double B, vector<int> Cobertura,  vector<double> &coef, vector<int> &clientes_corte, double &limitante, vector<double> &x, bool &corte_violado, double &y_til );
void LCIsequentiallifting(vector<double> Pesos , vector<int> Itens, double B, vector<int> Cobertura, vector<double> &coef, vector<int> &clientes_corte, double &limitante, vector<double> x, double &y_til );

bool verifica_corte_e_valido(vector<double> Pesos , double B, vector<double> coef, vector<int> clientes_corte, double limitante);

bool verifica_corte_e_violado(vector<double> Pesos , double B, vector<double> coef, vector<int> clientes_corte, double limitante, vector<double> x, double &y_til );

bool verifica_corte_e_violadoB(vector<double> Pesos , double B, vector<double> coef, vector<int> clientes_corte, double limitante, vector<double> x);

void LPIsuperaditivefunction (vector<double> Pesos , vector<int> Itens, double B, vector<int> vetor_pack, list<int> pack, vector<double> &alfa, vector<double> &lambda, vector<int> &vetor_clientes_restantes, double &limitante, vector<double> x, bool &corte_violado, double soma_demandas_pack );

void FCColumnGeneration(vector<double> Pesos , vector<int> Itens, double B, vector<double> &coef, vector<int> &clientes_corte, double &limitante, vector<double> &x, bool &corte_violado, bool &impossivel_cortar);

bool inteiroEstaNoVetor( vector<int>&vetor , int valor) ;

void ordena_clientes_crescente_demanda(vector<int>&vetor_clientes, vector<double> d);

bool vetorIdenticoNaMatriz( vector<vector<int>> matriz, vector<int> vetor, int k) ;

int maior_peso_possivel(const std::vector<int> items, int capacity);

int valor_minimo(const std::vector<int>& vetor1, const std::vector<int>& vetor2) ;

int indice_de_valor_minimo(const std::vector<int>& vetor1, const std::vector<int>& vetor2) ;

void imprime_vetor(vector<int> vetor) ;

void imprime_vetor_double( vector<double> vetor) ;

vector<int> cria_vetor_copia_de_lista(const std::list<int>& myList) ;

bool esta_na_lista(const std::list<int>& lista, int numero);

void cria_copia_lista (list<int> &lista_nova, list<int> lista);

void cria_copia_vetor_int(vector<int> &vetor_novo, vector<int> vetor);

void confere_solucao_s(int m, int n, vector<vector<double>> c,  vector<double> f, vector<list<int>>s, vector<double> d, list<int> abertas, list<int> fechadas, vector<double> pen, vector<double> p );

void confere_solucao_atende (int m, int n, vector<vector<double>> c,  vector<double> f, vector<int> atende, vector<double> d, vector<double> pen, double &fo_confere );

bool verifica_factivel_restrito(vector<int> vetor_facilidade, vector<double> pen );

void clientes_ordenados(int m, int n, vector<vector<double>> c, vector<vector<int>> &matriz_clientes_ordenados);

void facilidades_ordenadas(int m, int n, vector<vector<double>> c, vector<vector<int>> &matriz_facilidades_ordenadas);

void ordena_clientes_decrescente_demanda( vector<int>&vetor_clientes, vector<double> d);
void funcao_ordena_crescente(vector<int>&vetor_clientes, vector<double> d);
void funcao_ordena_decrescente(vector<int>&vetor_clientes, vector<double> d);
void ordena_clientes_crescente_demanda(vector<int>&vetor_clientes, vector<double> d);
void funcao_ordena_decrescente2(vector<int>&vetor_clientes, const vector<double>& x, const vector<double>& d);


void ordena_clientes_decrescente_demanda_2(vector<int>&vetor_clientes, vector<double> d);

bool verifica_factivel( vector<double> pen );

void imprime_lista(list<int> lista);

int verifica_lista_tabu_facilidades (vector<vector<int>> lista_tabu_facilidades, vector<int> vetor_facilidades_proximas1, int Q, int tamanho_lista_tabu_facilidades, int ultima_melhora);

void acrescenta_lista_tabu_facilidades (vector<vector<int>> &lista_tabu_facilidades, vector<int> cluster, int Q, int tamanho_lista_tabu_facilidades);

int verifica_lista_tabu_clientes (vector<vector<int>> lista_tabu_clientes, vector<int> vetor_conjunto_clientes, int Q, int tamanho_lista_tabu_clientes, int ultima_melhora);

void acrescenta_lista_tabu_clientes (vector<vector<int>> &lista_tabu_clientes, vector<int> vetor_conjunto_clientes, int tamanho_lista_clientes);

//void movimentos_fecha_abre_apos_solucao_inicial_best_improvement(double &fo, int m, int n, vector<vector<double>> c, vector<vector<int>> matriz_facilidades_ordenadas, vector<int> vetor_clientes, vector<double> p, vector<double> d, bool &conseguiu_alocar, list<int> &abertas_disponiveis, vector<double> &pen,  vector<int> &atende, vector<double> f, list<int> &fechadas_disponiveis, double menor_demanda, double soma_das_demandas, vector<int> &novas_facilidades,  vector<int> facilidades_potenciais);

void roda_relaxado_original (int m, int n, vector<double> f, vector<vector<double>> c, vector<double> d, vector<double> p, double soma_demandas, clock_t fim_CPU, clock_t inicio_CPU, double &tempo_relaxacao_no_raiz, double &limitante_inferior_no_raiz, int &quantidade_iteracoes, int &iteracao_que_encontrou_otimo, ofstream &arq_saida1, ofstream &arq_saida2, vector<list<int>>&corte1, vector<list<int>>&corte2, double &LB, double &primeiro_limitante_superior, double &otimo, vector<int> &atende, bool &factivel, vector<list<int>> &s, vector<double> &pen, list<int> &abertas, list<int> &fechadas, double &fo, list<int> &locais_promissores, vector<vector<int>>&N );

void roda_relaxado_original2 (int m, int n, vector<double> f, vector<vector<double>> c, vector<double> d, vector<double> p, double soma_demandas, clock_t fim_CPU, clock_t inicio_CPU, double &tempo_relaxacao_no_raiz, double &limitante_inferior_no_raiz, int &quantidade_iteracoes, int &iteracao_que_encontrou_otimo, ofstream &arq_saida1, ofstream &arq_saida2, vector<list<int>>&corte1, vector<list<int>>&corte2, double &LB, double &primeiro_limitante_superior, double &otimo, vector<int> &atende, bool &factivel, vector<list<int>> &s, vector<double> &pen, list<int> &abertas, list<int> &fechadas, double &fo, list<int> &locais_promissores,  vector<vector<int>>&N);

void verifica_se_problema_esparso_e_gap(list<int> locais_promissores, bool &problema_esparso_e_gap, vector<int> vetor_conjunto_clientes, vector<double> p, double soma_demandas, vector<int> vetor_facilidades_proximas1, double &soma_capacidades, vector<double> pen);

void CuttingPlane(int m, int n, vector<double> f, vector<double> p, vector<vector<double>> c, double &LB, double soma_demandas,   vector<double> d, vector<double> &y_copia, vector<vector<double>>  &x_copia, double capacidade_total, double demanda_total, clock_t &fim_CPU, clock_t &inicio_CPU , int &quantidade_iteracoes, int &iteracao_que_encontrou_otimo, double &tempo_resolve_problema_denso,  double tolerancia_pc,  int bigM, bool aplica_funcoes_de_lifting, bool aplica_lpi, bool aplica_lci_geral, bool aplica_fenchel_cuts,  list<vector<double>> &coeficientes_cortes_x, list<vector<int>> &clientes_cortes_x, vector<double> &limitantes_cortes_x, vector<int> &variavel_y_associada_ao_corte, list<vector<double>> &coeficientes_cortes_y, list<vector<int>> &facilidades_cortes_y, vector<double> &limitantes_cortes_y,  int max_iter_fc,  bool aplica_heuristica_gulosa1, bool aplica_heuristica_gulosa2, bool aplica_separacao_exata_de_cis, bool aplica_vubs, bool aplica_heuristica_construtiva, int tempo_maximo_pc, double tempo_limite_problema_estrategia_integralidade_parcial, IloEnv& env_denso, IloModel& mod_denso, IloCplex& cplex_denso, IloNumVarArray& y_denso, IloArray<IloNumVarArray> x_denso, bool &solucao_CPPIS_inteira, IloConstraintArray &restricoes_fcs );

bool vetorIgual(const std::vector<int>& a, const std::vector<int>& b);

void aplica_planos_cortantes_no_subproblema_cluster(vector<int> &vetor_conjunto_clientes, double &fo_antes,  double &fo_depois, vector<double> p, vector<vector<double>> c, vector<int> &atende_subsolucao, vector<double> d,  vector<int> &vetor_facilidades_proximas1, vector<double> f, bool coeficientes_inteiros, vector<int> atende, bool &cluster_gap, bool &melhorou, double fo,  int m, int n );

int verifica_lista_tabu_clientes (vector<vector<int>> lista_tabu_clientes, vector<int> vetor_conjunto_clientes, int Q, int tamanho_lista_tabu_clientes, int ultima_melhora);

int verifica_lista_tabu_novo(vector<vector<int>> lista_tabu, vector<int> vetor, int inicio_da_lista, int final_da_lista); 

void acrescenta_lista_tabu_clientes (vector<vector<int>> &lista_tabu_clientes, vector<int> vetor_conjunto_clientes, int tamanho_lista_clientes);




#define MAXSTATES 400000 

#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <stdarg.h>
#include <values.h>
#include <string.h>
#include <math.h>
#include <malloc.h>
#include <cstdlib>


/* ======================================================================
				   macros
   ====================================================================== */

#define SYNC            5      /* when to switch to linear scan in bins */
#define SORTSTACK     200      /* depth of stack used in qsort */
#define MINMED        100      /* find exact median in qsort if larger size */

#define TRUE  1
#define FALSE 0

#define LEFT  1
#define RIGHT 2

#define PARTIATE 1
#define SORTALL  2

#define MAXV (8*sizeof(btype)) /* number of bits in a long integer */
#define PMAX 1                 /* profit of worlds most efficient item  */
#define WMAX 0                 /* weight of worlds most efficient item  */
#define PMIN 0                 /* profit of worlds least efficient item */
#define WMIN 1                 /* weight of worlds least efficient item */

#define DET(a1, a2, b1, b2)    ((a1) * (ptype) (b2) - (a2) * (ptype) (b1))
#define SWAP(a, b)   { register item t; t = *(a); *(a) = *(b); *(b) = t; }
#define DIFF(a,b)              ((int) ((b)-(a)+1))
#define NO(a,p)                ((int) ((p) - (a)->fitem + 1))
#define N(a,p)                 ((int) ((p) - (a)->d.set1))
#define L(x)                   ((long) (x))
#define SZ(a)                  (*(((int *) (a)) - 4) - 1)


/* ======================================================================
				 type declarations
   ====================================================================== */

typedef int           boolean;
typedef long          ntype;   /* number of states/items   */
typedef long          itype;   /* item profits and weights */
typedef long          stype;   /* sum of pofit or weight   */
typedef double        ptype;   /* product type (sufficient precision) */
typedef unsigned long btype;   /* binary representation of solution */

/* item record */
typedef struct irec {
  itype   p;     /* profit */
  itype   w;     /* weight */
  boolean *x;    /* solution variable */
} item;

typedef struct { /* i-stack */
  item  *f;      /* first item in interval */
  item  *l;      /* last item in interval */
} interval;

/* state in dynamic programming */
typedef struct pv {
  stype psum;    /* profit sum */
  stype wsum;    /* weight sum */
  btype vect;    /* solution vector */
} state;

/* set of states */
typedef struct pset {
  ntype size;    /* set size */
  state *fset;   /* first element in set */
  state *lset;   /* last element in set */
  state *set1;   /* first element in array */
  state *setm;   /* last element in array */
} stateset;

typedef struct { /* all problem information */
  ntype    n;               /* number of items         */
  item     *fitem;          /* first item in problem   */
  item     *litem;          /* last item in problem    */
  item     *ftouch;         /* first item considered for reduction */
  item     *ltouch;         /* last item considered for reduction */
  item     *s;              /* current core is [s,t]   */
  item     *t;              /*                         */
  item     *b;              /* break item              */
  item     *fpart;          /* first item returned by partial sort */
  item     *lpart;          /* last item returned by partial sort */
  stype    wfpart;          /* weight sum up to fpart  */
  item     *fsort;          /* first sorted item       */
  item     *lsort;          /* last sorted item        */
  stype    wfsort;          /* weight sum up to fsort  */
  stype    c;               /* current capacity        */
  stype    cstar;           /* origianl capacity       */
  stype    z;               /* current solution        */
  stype    zstar;           /* optimal solution        */
  stype    zwsum;           /* weight sum of zstar     */
  itype    ps, ws, pt, wt;  /* items for deriving bounds */

  btype    vno;             /* current vector number   */
  item *   vitem[MAXV];     /* current last MAXV items */
  item *   ovitem[MAXV];    /* optimal set of items    */
  btype    ovect;           /* optimal solution vector */

  stype    dantzig;         /* dantzig upper bound     */
  stype    ub;              /* global upper bound      */
  stype    psumb;           /* profit sum up to b      */
  stype    wsumb;           /* weight sum up to b      */
  boolean  firsttime;       /* used for restoring x    */
  boolean  welldef;         /* is x welldefined        */
  stateset  d;              /* set of partial vectors  */
  interval *intv1, *intv2;
  interval *intv1b, *intv2b;

  /* debug */
  long     iterates;        /* counters used to obtain specific */
  long     simpreduced;     /* information about the solution process */
  long     pireduced;
  long     pitested;
  long     maxstates;
  long     coresize;
  long     bzcore;
} allinfo;


void errorx(char *str, ...);
void pfree(void *p);
void *palloc(long size);
state *findvect(stype ws, state *f, state *l);
void push(allinfo *a, int side, item *f, item *l);
void improvesolution(allinfo *a, state *v);
void definesolution(allinfo *a);
item *median(item *f1, item *l1, ntype s);
void partsort(allinfo *a, item *f, item *l, stype ws, int what);
boolean haschance(allinfo *a, item *i, int side);
void multiply(allinfo *a, item *h, int side);
void simpreduce(int side, item **f, item **l, allinfo *a);
void reduceset(allinfo *a);
void initfirst(allinfo *a, stype ps, stype ws);
void initvect(allinfo *a);
void copyproblem(item *f, item *l, int *p, int *w, int *x);
void findbreak(allinfo *a);
stype minknap(int n, int *p, int *w, int *x, int c);
