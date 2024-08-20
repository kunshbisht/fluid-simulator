#include "raylib.h"
#include "raymath.h"
#include <iostream>

const int N = 128;
const int iter = 10;
const int scale = 4;
const int fps = 500;

// Updated IX function for 3D grids
int IX(int x, int y) {
    return x + y * N;
}

class Fluid {
    public:
        int size = N;
        float dt;
        float diff;
        float visc;

        float s[N*N];
        float density[N*N];

        float Vx[N*N];
        float Vy[N*N];
        float Vx0[N*N];
        float Vy0[N*N];

        void step() {
            diffuse(1, Vx0, Vx, visc);
            diffuse(2, Vy0, Vy, visc);
            
            project(Vx0, Vy0, Vx, Vy);
            
            advect(1, Vx, Vx0, Vx0, Vy0);
            advect(2, Vy, Vy0, Vx0, Vy0);
            
            project(Vx, Vy, Vx0, Vy0);
            
            diffuse(0, s, density, diff);
            advect(0, density, s, Vx, Vy);
        }

        void addDensity(int x, int y, float amount) {
            int index = IX(x, y);
            density[index] += amount;
        }

        void addVelocity(int x, int y, float amountX, float amountY) {
            int index = IX(x, y);
            Vx[index] += amountX;
            Vy[index] += amountY;
        }

        void diffuse(int b, float x[N*N], float x0[N*N], float diff) {
            float a = dt * diff * (N-2) * (N-2);
            lin_solve(b, x, x0, a, 1 + 6 * a);
        }

        void set_bnd(int b, float x[]) {
            for(int i = 1; i < N - 1; i++) {
                x[IX(i, 0  )] = b == 2 ? -x[IX(i, 1  )] : x[IX(i, 1  )];
                x[IX(i, N-1)] = b == 2 ? -x[IX(i, N-2)] : x[IX(i, N-2)];
            }

            for(int j = 1; j < N - 1; j++) {
                x[IX(0  , j)] = b == 1 ? -x[IX(1  , j)] : x[IX(1  , j)];
                x[IX(N-1, j)] = b == 1 ? -x[IX(N-2, j)] : x[IX(N-2, j)];
            }

            x[IX(0, 0)]     = 0.5f * (x[IX(1, 0)]
                                + x[IX(0, 1)]);
            x[IX(N-1, 0)]   = 0.5f * (x[IX(N-2, 0)]
                                + x[IX(N-1, 1)]);
            x[IX(0, N-1)]   = 0.5f * (x[IX(1, N-1)]
                                + x[IX(0, N-2)]);
            x[IX(N-1, N-1)] = 0.5f * (x[IX(N-2, N-1)]
                                + x[IX(N-1, N-2)]);
        }

        void lin_solve(int b, float x[], float x0[], float a, float c) {
            float cRecip = 1.0 / c;
            for (int k = 0; k < iter; k++) {
                for (int j = 1; j < N - 1; j++) {
                    for (int i = 1; i < N - 1; i++) {
                        x[IX(i, j)] = (x0[IX(i, j)] + a * (x[IX(i+1, j)] + x[IX(i-1, j)] + x[IX(i, j+1)] + x[IX(i, j-1)])) * cRecip;
                    }
                }
                set_bnd(b, x);
            }
        }

        void project(float velocX[], float velocY[], float p[], float div[]) {
            for (int j = 1; j < N - 1; j++) {
                for (int i = 1; i < N - 1; i++) {
                    div[IX(i, j)] = -0.5f * (velocX[IX(i+1, j)] - velocX[IX(i-1, j)] + velocY[IX(i, j+1)] - velocY[IX(i, j-1)]) / N;
                    p[IX(i, j)] = 0;
                }
            }
            set_bnd(0, div); 
            set_bnd(0, p);
            lin_solve(0, p, div, 1, 6);
            
            for (int j = 1; j < N - 1; j++) {
                for (int i = 1; i < N - 1; i++) {
                    velocX[IX(i, j)] -= 0.5f * (p[IX(i+1, j)] - p[IX(i-1, j)]) * N;
                    velocY[IX(i, j)] -= 0.5f * (p[IX(i, j+1)] - p[IX(i, j-1)]) * N;
                }
            }
            set_bnd(1, velocX);
            set_bnd(2, velocY);
        }

        void advect(int b, float d[], float d0[],  float velocX[], float velocY[]) {
            float i0, i1, j0, j1;
            
            float dtx = dt * (N - 2);
            float dty = dt * (N - 2);
            
            float s0, s1, t0, t1, u0, u1;
            float tmp1, tmp2, x, y;
            
            float Nfloat = N;
            float ifloat, jfloat;
            int i, j;
            
            for(j = 1, jfloat = 1; j < N - 1; j++, jfloat++) { 
                for(i = 1, ifloat = 1; i < N - 1; i++, ifloat++) {
                    tmp1 = dtx * velocX[IX(i, j)];
                    tmp2 = dty * velocY[IX(i, j)];
                    x    = ifloat - tmp1; 
                    y    = jfloat - tmp2;
                    
                    if(x < 0.5f) x = 0.5f; 
                    if(x > Nfloat + 0.5f) x = Nfloat + 0.5f; 
                    i0 = floorf(x); 
                    i1 = i0 + 1.0f;
                    if(y < 0.5f) y = 0.5f; 
                    if(y > Nfloat + 0.5f) y = Nfloat + 0.5f; 
                    j0 = floorf(y);
                    j1 = j0 + 1.0f; 
                    
                    s1 = x - i0; 
                    s0 = 1.0f - s1; 
                    t1 = y - j0; 
                    t0 = 1.0f - t1;
                    
                    int i0i = i0;
                    int i1i = i1;
                    int j0i = j0;
                    int j1i = j1;
                    
                    d[IX(i, j)] = 
                        s0 * ( t0 * d0[IX(i0i, j0i)] + t1 * d0[IX(i0i, j1i)])
                        + s1 * ( t0 * d0[IX(i1i, j0i)] + t1 * d0[IX(i1i, j1i)]);
                }
            }
            set_bnd(b, d);
        }

        void renderD() {
            for (int i = 0; i < N; i++) {
                for (int j = 0; j < N; j++) {
                    float x = i * scale;
                    float y = j * scale;
                    float d = density[IX(i, j)];

                    // Map density to an intensity value for hue (0 to 360 degrees)
                    float hue = d;
                    float saturation = 1.0f; // Full saturation
                    float value = 1.0f;      // Full brightness

                    // Convert HSV to RGB
                    Color color = ColorFromHSV(hue, saturation, value);

                    // Draw rectangle with the HSV-based color
                    DrawRectangle(x, y, scale, scale, color);
                }
            }
        }
};

Fluid fluid;

int main() {
    InitWindow(N*scale, N*scale, "Raylib Fluid Simulation");

    SetTargetFPS(fps);

    fluid = Fluid();
    fluid.dt = 1.0f/fps;
    fluid.diff = 0;
    fluid.visc = 0;

    // Main game loop
    while (!WindowShouldClose()) {
        // Update logic here

        // Draw
        BeginDrawing();
        ClearBackground(BLACK);

        fluid.step();

        if(IsMouseButtonDown(0)) {
            if (std::max(GetMouseX(), GetMouseY()) <= N * scale && std::min(GetMouseX(), GetMouseY()) >= 0) {
                fluid.addDensity(GetMouseX()/scale, GetMouseY()/scale, 500);
                Vector2 vel = Vector2Scale(Vector2Normalize(GetMouseDelta()), 10);
                fluid.addVelocity(GetMouseX()/scale, GetMouseY()/scale, vel.x, vel.y);
            }
        }

        fluid.renderD();

        EndDrawing();
    }

    // De-Initialization
    CloseWindow(); 

    return 0;
}