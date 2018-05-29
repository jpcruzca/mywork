#include <vector>

//constants
const double u = 0.012277471;
const double T = 17.0652165601579625588917206249;

//function declaration
void initial_condition(std::vector<double> & pos, std::vector<double> & vel);

void print(const std::vector<double> & v);

void rk4(std::vector<double> & pos, std::vector<double> & vel, const double tini, const double tend);

double compute(const std::vector<double> & pos, const std::vector<double> & vel, const double t, const int id); // retorna valores para la posicion, velocidad, aceleracion

double rvar(const std::vector<double> & pos, const std::vector<double> & vel);

double svar(const std::vector<double> & pos, const std::vector<double> & vel);
void init_gnuplot(void);
void plot_gnuplot(const std::vector<double>& m);
