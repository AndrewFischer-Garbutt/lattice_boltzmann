#include <vector>
#include <iostream>
#include <algorithm>
#include <iomanip>      // std::setprecision
#include "graphics.h"

#define GRID_WIDTH 100
#define GRID_HIGHT 100
#define NUM_COMPONENTS 9
#define tau 1.1

class LBM
{
	private:
		struct I2
		{
			I2(int x_in, int y_in)
				: x(x_in), y(y_in) {}

			I2()
				: x(0), y(0) {}

			int x;
			int y;
		};

		struct F2
		{
			F2(float x_in, float y_in)
				: x(x_in), y(y_in) {}

			F2()
				: x(0), y(0) {}

			float x;
			float y;
		};

		float f[GRID_HIGHT][GRID_WIDTH][NUM_COMPONENTS];
		float f_star[GRID_HIGHT][GRID_WIDTH][NUM_COMPONENTS];
		float f_eq[GRID_HIGHT][GRID_WIDTH][NUM_COMPONENTS];
		float rho[GRID_HIGHT*GRID_WIDTH];
		F2 u[GRID_HIGHT][GRID_WIDTH];
		I2 e[9] = { {0,0},{1,0},{0,-1},{-1,0},{0,1},{1,-1},{-1,-1},{-1,1},{1,1} };
		float w[9] = { 4./9, 1./9, 1./9, 1./9, 1./9, 1./36, 1./36, 1./36, 1./36 };
		float c = 1.;

	public:
		// Default construct to shut up compiler errors.
		LBM() {}

		void initialize()
		{
			for(int i = 0; i < GRID_HIGHT; i++)
			{
				for(int j = 0; j < GRID_WIDTH; j++)
				{
					rho[i+GRID_HIGHT*j] = 0.;
					for(int k = 0; k < NUM_COMPONENTS; k++)
					{
						if(25 <= i && i <= 50 && 25 <= j && j <= 50 && k == 1) f[i][j][k] = 2./10;
						else if(25 <= i && i <= 50 && 25 <= j && j <= 50) f[i][j][k] = 1./10;
						else f[i][j][k] = 1./9;
					}
				}
			}
		}

		void streaming_step()
		{
			for(int i = 0; i < GRID_HIGHT; i++)
			{
				for(int j = 0; j < GRID_WIDTH; j++)
				{
					for(int k = 0; k < NUM_COMPONENTS; k++)
					{
						int x_pos = ((j + e[k].x)%GRID_WIDTH + GRID_WIDTH) % GRID_WIDTH;
						int y_pos = ((i + e[k].y)%GRID_HIGHT + GRID_HIGHT) % GRID_HIGHT;
						f_star[y_pos][x_pos][k] = f[i][j][k];
					}
				}
			}
		}

		void compute_macroscopic()
		{
			for(int i = 0; i < GRID_HIGHT; i++)
			{
				for(int j = 0; j < GRID_WIDTH; j++)
				{
					rho[i+GRID_HIGHT*j] = 0.;
					for(int k = 0; k < NUM_COMPONENTS; k++)
					{
						rho[i+GRID_HIGHT*j] += f_star[i][j][k];
					}
					u[i][j] = {0.,0.};
					for(int k = 0; k < NUM_COMPONENTS; k++)
					{
						u[i][j].x += c*f_star[i][j][k]*e[k].x / rho[i+GRID_HIGHT*j];
						u[i][j].y += c*f_star[i][j][k]*e[k].y / rho[i+GRID_HIGHT*j];
					}
				}
			}
		}

		void update_f_eq()
		{
			for(int i = 0; i < GRID_HIGHT;i++)
			{
				for(int j = 0; j < GRID_WIDTH; j++)
				{
					float u_mag2 = u[i][j].x * u[i][j].x + u[i][j].y * u[i][j].y;
					for(int k = 0; k < NUM_COMPONENTS; k++)
					{
						float dot_p = e[k].x * u[i][j].x + e[k].y * u[i][j].y;

						float s = w[k]*(3.*dot_p/c + (9./2)*(dot_p*dot_p)/(c*c) - (3./2)*(u_mag2)/(c*c));
						f_eq[i][j][k] = rho[i+GRID_HIGHT*j]*(s + w[k]);
					}
				}
			}
		}

		void update_f()
		{
			for(int i = 0; i < GRID_HIGHT; i++)
				for(int j = 0; j < GRID_WIDTH; j++)
					for(int k = 0; k < NUM_COMPONENTS; k++)
						f[i][j][k] = f_star[i][j][k] - (1. / tau)*(f_star[i][j][k] - f_eq[i][j][k]);
		}

		void set_up_display()
		{
			Graphics::setup(800, 800, GRID_HIGHT, GRID_WIDTH);
		}

		void update_display()
		{
			float mx = -1, mn = std::numeric_limits<int>::max();
			for(int y = 0; y < GRID_HIGHT; ++y)
				for(int x = 0; x < GRID_WIDTH; ++x)
				{
					if(rho[x+y*GRID_HIGHT] > mx) mx = rho[y+x*GRID_HIGHT];
					if(rho[x+y*GRID_HIGHT] < mn) mn = rho[y+x*GRID_HIGHT];
				}

			float color[GRID_HIGHT*GRID_WIDTH];
			for(int y = 0; y < GRID_HIGHT; ++y)
				for(int x = 0; x < GRID_WIDTH; ++x)
					color[(GRID_HIGHT - 1 - y) + GRID_HIGHT*x] = (rho[y+x*GRID_HIGHT] - mn) / (mx - mn);
			Graphics::draw(color);
		}
};

int main()
{
	LBM* lbm = new LBM();
	lbm->initialize();
	lbm->set_up_display();

	while(true)
	{
		lbm->update_display();
		lbm->streaming_step();
		lbm->compute_macroscopic();
		lbm->update_f_eq();
		lbm->update_f();
	}

	delete lbm;
}
