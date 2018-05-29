/*#include <vector>
//constantes
const double u=0.012277471; //proporcion entre las masas
const double T=17.0652165601579625588917206249; //periodo
const double eps=0.001;

//functions
void initial_condition(std::vector<double> & pos, std::vector<double> & vel);

void print(const std::vector<double> & v);

void dp54(std::vector<double> & pos, std::vector<double> & vel, const double tini, const double tend);

double compute(const std::vector<double> & pos, const std::vector<double> & vel, const double t, const int id); // retorna valores de posicion, velocidad, aceleracion

double rvar(const std::vector<double> & pos, const std::vector<double> & vel);

double svar(const std::vector<double> & pos, const std::vector<double> & vel);

//automatic step-size control functions


double saux(const double x, const double dt);

double norm(const std::vector<double> & x);
*/

#include <vector>
const double u=0.012277471; //proporcion entre las masas
const double T=17.0652165601579625588917206249; //periodo
const double eps=0.001;
//functions
void initial_condition(std::vector<double> & pos, std::vector<double> & vel);

void print(const std::vector<double> & v);

void dp54(std::vector<double> & pos, std::vector<double> & vel, const double tini, const double tend);

double compute(const std::vector<double> & pos, const std::vector<double> & vel, const double t, const int id); // retorna valores de posicion, velocidad, aceleracion

double rvar(const std::vector<double> & pos, const std::vector<double> & vel);

double svar(const std::vector<double> & pos, const std::vector<double> & vel);

//automatic step-size control functions
//double dtnew(double e,double dt);
double dtnew(const double p, const double v, const double dt);
double saux(const double x, const double dt);

double norm(const std::vector<double> & x);

double E(const std::vector<double> & vel, const std::vector<double> & velaux,const std::vector<double> & pos, const std::vector<double> & posaux);//medida conjunta del error

