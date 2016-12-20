#include "Sde.h"

Sde::Sde(string solution_file)
:Ode(solution_file)
{
    diagonal=false;
    commutative=false;
}

Sde::~Sde()
{
}

void Sde::reinit()
{
    init_solution();
    t=t0;
}

int Sde::brownian_size()
{
    return Wsize;
}

bool Sde::is_commutative() const
{
    return commutative;
}

bool Sde::isDiagonal() const
{
    return diagonal;
}
void Sde::sample(Real t, Real h)
{
    W->sample(t,h);
}

void Sde::setW(BrownianMotion* W)
{
    this->W = W;
    Wsize = W->size();
}

Vector& Sde::getIr()
{
    return W->getIr();
}

Vector& Sde::getChir()
{
    return W->getChir();
}

Vector& Sde::getIrr()
{
    return W->getIrr();
}

vector<Vector>& Sde::getIpq()
{
    return W->getIpq();
}

Vector& Sde::getOnes()
{
    return W->getOnes();
}

void Sde::g(Real t, Vector& x, Vector& G, int r)
{
    cout<<"ERROR: using undefined function g(Real, Vector&, Vector&, int)"<<endl;
}

void Sde::g(Real t, Vector& x, Vector& G)
{
    cout<<"ERROR: using undefined function g(Real, Vector&, Vector&)"<<endl;
}

void Sde::g(Real t, Vector& x, Matrix& G)
{
    cout<<"ERROR: using undefined function g(Real, Vector&, Matrix&)"<<endl;
}