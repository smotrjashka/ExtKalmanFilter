// KalmanFilter2.cpp : This file contains the 'main' function. Program execution begins and ends there.
//
#define SDL_MAIN_HANDLED
#include <iostream>
#include <cmath>
#include <limits>
#include <SDL2/SDL.h>
#include <SDL2/SDL_ttf.h>
#include <algorithm>
#include <thread>


long double func_fx(long double mu, long double sigma2, long double x);

long double func_new_mean(long double prior_mean, long double prior_sigma2, long double measured_mean, long double measured_sigma2);

long double func_new_sigma2(long double prior_sigma2, long double measured_sigma2);


using namespace std;


const double R_STANDART_DEVIATION = 0.1;  //measurement sigma
const double PSI_STANDART_DEVIATION = 3;
SDL_Renderer* renderer;


struct gausian {
    long double mean;
    long double variance;

    gausian(long double mean, long double variance) : mean(mean), variance(variance) {}
};

class State_Vector {
public:
    long double x;
    long double y;
    long double vx;
    long double vy;

    State_Vector(long double x, long double y, long double vx, long double vy) : x(x), y(y), vx(vx), vy(vy) {}

    State_Vector(long double x, long double y) : x(x), y(y) { vx = 0; vy = 0; }

    State_Vector() { x = 0; y = 0; vx = 0; vy = 0; }

};

template <int N, int M>
class Matrix {
public:
    double matr[N][M]{};

    void print_string() {
        for (int i = 0; i < N; ++i) {
            for (int j = 0; j < M; ++j) {
                double elem = matr[i][j];
                cout << elem << ", ";
            }
            cout << endl;
        }
        cout << endl;
    }
};

template <int N>
struct State {
    State() {}

    State(Matrix<4, 1> matrix, Matrix<4, 4> matrix1) {
        state_vector = matrix;
        covariance_matrix = matrix1;
    }

    Matrix<N, 1> state_vector;
    Matrix<N, N> covariance_matrix;
};

template <int N, int M, int K>
Matrix<N, K> multipleMatrix(Matrix<N, M> m1, Matrix<M, K> m2) {
    Matrix<N, K> resultMatrix;

    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < K; ++j) {
            for (int k = 0; k < M; ++k) {
                resultMatrix.matr[i][j] += m1.matr[i][k] * m2.matr[k][j];
            }

        }
    }

    return resultMatrix;
}

template <int N, int M>
Matrix<N, M> addMatrix(Matrix<N, M> m1, Matrix<N, M> m2) {
    Matrix<N, M> resultMatrix;

    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < M; ++j) {

            resultMatrix.matr[i][j] = m1.matr[i][j] + m2.matr[i][j];
        }
    }

    return resultMatrix;
}

template <int N, int M>
Matrix<N, M> minusMatrix(Matrix<N, M> m1, Matrix<N, M> m2) {
    Matrix<N, M> resultMatrix;

    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < M; ++j) {

            resultMatrix.matr[i][j] = m1.matr[i][j] - m2.matr[i][j];
        }
    }

    return resultMatrix;
}

template <int N, int M>
Matrix<N, M> addMatrix(double a, Matrix<N, M> m1, double b, Matrix<N, M> m2) {
    Matrix<N, M> resultMatrix;

    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < M; ++j) {

            resultMatrix.matr[i][j] = a * m1.matr[i][j] + b * m2.matr[i][j];
        }
    }

    return resultMatrix;
}

template <int N, int M>
Matrix<N, M> multipleMatrix(double a, Matrix<N, M> m1) {
    Matrix<N, M> resultMatrix;

    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < M; ++j) {

            resultMatrix.matr[i][j] = a * m1.matr[i][j];
        }
    }

    return resultMatrix;
}

template <int N, int M>
Matrix<M, N> transposeMatrix(Matrix<N, M> m1) {
    Matrix<M, N> resultMatrix;

    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < M; ++j) {

            resultMatrix.matr[j][i] = m1.matr[i][j];
        }
    }

    return resultMatrix;
}

Matrix<4, 4> get_motion_matrix(double t) {
    Matrix<4, 4> result;

    for (int i = 0; i < 4; ++i) {
        result.matr[i][i] = 1;
    }

    result.matr[0][2] = t;
    result.matr[1][3] = t;

    return result;
}

Matrix<4, 4> get_velocity_noise_covariance_matrix(double q, double t) {
    Matrix<4, 4> result;

    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 4 - i; ++j) {
            if ((i + j) % 2 != 0) {
                result.matr[i][j] = 0;
                result.matr[j][i] = 0;
            }
            else if (abs(i - j) == 2) {

                result.matr[i][j] = pow(t, 2) / 2;
                result.matr[j][i] = result.matr[i][j];
            }
        }
    }

    result.matr[0][0] = pow(t, 3) / 3;
    result.matr[1][1] = pow(t, 3) / 3;
    result.matr[2][2] = t;
    result.matr[3][3] = t;

    if (q != 1) {
        return multipleMatrix(q, result);
    }

    return result;
}

Matrix<4, 4> get_diff_t_velocity_covariance_matrix(double q, double t) {
    Matrix<4, 4> differByT;

    double mas[4][4] = { {pow(t, 2), 0, t, 0}, {0, pow(t, 2), 0, t}, {t, 0, 1, 0}, {0, t, 0, 1} };
    differByT = reinterpret_cast<Matrix<4, 4>&&>(mas);

    if (q != 1) {
        return multipleMatrix(q, differByT);
    }

    return differByT;
}

//we asume that noise_mean = 0
Matrix<4, 1> prediction_step_state_vector(Matrix<4, 1> initial_state_vector,/* Matrix<4, 1> noise_mean,*/ Matrix<4, 4> motion_matrix) {

    Matrix<4, 1> result_matrix;

    result_matrix = multipleMatrix(motion_matrix, initial_state_vector);
    //   result_matrix = addMatrix(result_matrix, noise_mean);

    return result_matrix;

}

//motion_matrix for Extended is real jacob (derivative) of motion matrix
Matrix<4, 4> predict_step_covar_matrix(Matrix<4, 4> motion_matrix, Matrix<4, 4> prev_predict_covar, Matrix<4, 4> Q_covar) {
    Matrix<4, 4> result;

    result = addMatrix(multipleMatrix(multipleMatrix(motion_matrix, prev_predict_covar), transposeMatrix(motion_matrix)),
        Q_covar);

    return result;
}


Matrix<4, 4> getMeasurmentModelMatrix(/*double dt*/) {
    Matrix<4, 4> measur_model;

    measur_model.matr[0][0] = 1;
    measur_model.matr[1][1] = 1;
    measur_model.matr[0][2] = 0; //dt;
    measur_model.matr[1][3] = 0; /*dt;*/

    return measur_model;
}

/*Matrix<2, 4> get_diff_measured_model_matrix(double dt){
    Matrix<2, 4> meas_model;

    meas_model.matr[0][2] = 1;
    meas_model.matr[1][3] = 1;

    return meas_model;
}*/

/*Matrix<4, 4> deltaH_numeric(Matrix<4, 1> predicted_state, Matrix<2, 1> measured_state){
    Matrix<4, 4> delta_h;

    return delta_h;
}*/


double matrixDet(Matrix<3, 3> matrix) {
    double determ = 0.0;

    determ = matrix.matr[0][0] * matrix.matr[1][1] * matrix.matr[2][2]
        + matrix.matr[0][1] * matrix.matr[1][2] * matrix.matr[2][0]
        + matrix.matr[0][2] * matrix.matr[1][0] * matrix.matr[2][1]
        - matrix.matr[0][2] * matrix.matr[1][1] * matrix.matr[2][0]
        - matrix.matr[0][0] * matrix.matr[1][2] * matrix.matr[2][1]
        - matrix.matr[0][1] * matrix.matr[1][0] * matrix.matr[2][2];

    return determ;
}

double matrixDet(Matrix<4, 4> matrix) {
    double determ = 0.0;
    Matrix<3, 3> minor_matrix0, minor_matrix1, minor_matrix2, minor_matrix3;

    for (int i = 1; i < 4; i++)
    {
        int j0 = 0, j1 = 0, j2 = 0, j3 = 0;
        for (int j = 0; j < 4; j++)
        {
            double temp_value = matrix.matr[i][j];
            if (j != 0) {
                minor_matrix0.matr[i - 1][j0] = temp_value;
                j0++;
            }
            if (j != 1) {
                minor_matrix1.matr[i - 1][j1] = temp_value;
                j1++;
            }

            if (j != 2) {
                minor_matrix2.matr[i - 1][j2] = temp_value;
                j2++;
            }

            if (j != 3) {
                minor_matrix3.matr[i - 1][j3] = temp_value;
                j3++;
            }
        }
    }

    determ = matrix.matr[0][0] * matrixDet(minor_matrix0) - matrix.matr[0][1] * matrixDet(minor_matrix1)
        + matrix.matr[0][2] * matrixDet(minor_matrix2) - matrix.matr[0][3] * matrixDet(minor_matrix3);

    return determ;
}


double matrixDet(Matrix<2, 2> matrix) {
    double determ = 0.0;

    determ = matrix.matr[0][0] * matrix.matr[1][1] - matrix.matr[0][1] * matrix.matr[1][0];

    return determ;
}

Matrix<4, 4> reverse_matrix(Matrix<4, 4> matrix_for_reverse) {
    Matrix<4, 4> result;

    //  cout << "tart reverse this matrix: " << endl;
   //   matrix_for_reverse.print_string();
    double determinant = 0.0;
    determinant = matrixDet(matrix_for_reverse);


    Matrix<3, 3> matrix00;
    for (int i = 1; i < 4; ++i) {
        for (int j = 1; j < 4; ++j) {
            matrix00.matr[i - 1][j - 1] = matrix_for_reverse.matr[i][j];
        }
    }
    // matrix00.print_string();
    result.matr[0][0] = matrixDet(matrix00);

    Matrix<3, 3> matrix01;
    for (int i = 1; i < 4; ++i) {
        for (int j = 0; j < 4; ++j) {
            if (j < 1) {
                matrix01.matr[i - 1][j] = matrix_for_reverse.matr[i][j];
            }
            else if (j > 1) {
                matrix01.matr[i - 1][j - 1] = matrix_for_reverse.matr[i][j];
            }
        }
    }
    // matrix01.print_string();
    result.matr[0][1] = -1 * matrixDet(matrix01);

    Matrix<3, 3> matrix02;
    for (int i = 1; i < 4; ++i) {
        for (int j = 0; j < 4; ++j) {
            if (j < 2) {
                matrix02.matr[i - 1][j] = matrix_for_reverse.matr[i][j];
            }
            else if (j > 2) {
                matrix02.matr[i - 1][j - 1] = matrix_for_reverse.matr[i][j];
            }
        }
    }
    // matrix02.print_string();
    result.matr[0][2] = matrixDet(matrix02);

    Matrix<3, 3> matrix03;
    for (int i = 1; i < 4; ++i) {
        for (int j = 0; j < 3; ++j) {
            matrix03.matr[i - 1][j] = matrix_for_reverse.matr[i][j];
        }
    }
    // matrix03.print_string();
    result.matr[0][3] = -1 * matrixDet(matrix03);

    Matrix<3, 3> matrix10;
    for (int i = 0; i < 4; ++i) {

        for (int j = 1; j < 4; ++j) {
            if (i < 1) {
                matrix10.matr[i][j - 1] = matrix_for_reverse.matr[i][j];
            }
            else if (i > 1) {
                matrix10.matr[i - 1][j - 1] = matrix_for_reverse.matr[i][j];
            }
        }
    }
    //   matrix10.print_string();
    result.matr[1][0] = (-1) * matrixDet(matrix10);

    Matrix<3, 3> matrix11;
    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 4; ++j) {
            if (j < 1) {
                if (i < 1) {
                    matrix11.matr[i][j] = matrix_for_reverse.matr[i][j];
                }
                else if (i > 1) {
                    matrix11.matr[i - 1][j] = matrix_for_reverse.matr[i][j];
                }
            }
            else if (j > 1) {
                if (i < 1) {
                    matrix11.matr[i][j - 1] = matrix_for_reverse.matr[i][j];
                }
                else if (i > 1) {
                    matrix11.matr[i - 1][j - 1] = matrix_for_reverse.matr[i][j];
                }
            }
        }
    }
    //  matrix11.print_string();
    result.matr[1][1] = matrixDet(matrix11);


    Matrix<3, 3> matrix12;
    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 4; ++j) {
            if (j < 2) {
                if (i < 1) {
                    matrix12.matr[i][j] = matrix_for_reverse.matr[i][j];
                }
                else if (i > 1) {
                    matrix12.matr[i - 1][j] = matrix_for_reverse.matr[i][j];
                }
            }
            else if (j > 2) {
                if (i < 1) {
                    matrix12.matr[i][j - 1] = matrix_for_reverse.matr[i][j];
                }
                else if (i > 1) {
                    matrix12.matr[i - 1][j - 1] = matrix_for_reverse.matr[i][j];
                }
            }
        }
    }
    //   matrix12.print_string();
    result.matr[1][2] = -1 * matrixDet(matrix12);

    Matrix<3, 3> matrix13;
    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 3; ++j) {
            if (i < 1) {
                matrix13.matr[i][j] = matrix_for_reverse.matr[i][j];
            }
            else if (i > 1) {
                matrix13.matr[i - 1][j] = matrix_for_reverse.matr[i][j];
            }
        }
    }
    //   matrix13.print_string();
    result.matr[1][3] = matrixDet(matrix13);

    Matrix<3, 3> matrix20;
    for (int i = 0; i < 4; ++i) {

        for (int j = 1; j < 4; ++j) {
            if (i < 2) {
                matrix20.matr[i][j - 1] = matrix_for_reverse.matr[i][j];
            }
            else if (i > 2) {
                matrix20.matr[i - 1][j - 1] = matrix_for_reverse.matr[i][j];
            }
        }
    }
    // matrix20.print_string();
    result.matr[2][0] = matrixDet(matrix20);

    Matrix<3, 3> matrix21;
    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 4; ++j) {
            if (j < 1) {
                if (i < 2) {
                    matrix21.matr[i][j] = matrix_for_reverse.matr[i][j];
                }
                else if (i > 2) {
                    matrix21.matr[i - 1][j] = matrix_for_reverse.matr[i][j];
                }
            }
            else if (j > 1) {
                if (i < 2) {
                    matrix21.matr[i][j - 1] = matrix_for_reverse.matr[i][j];
                }
                else if (i > 2) {
                    matrix21.matr[i - 1][j - 1] = matrix_for_reverse.matr[i][j];
                }
            }
        }
    }
    //   matrix21.print_string();
    result.matr[2][1] = -1 * matrixDet(matrix21);


    Matrix<3, 3> matrix22;
    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 4; ++j) {
            if (j < 2) {
                if (i < 2) {
                    matrix22.matr[i][j] = matrix_for_reverse.matr[i][j];
                }
                else if (i > 2) {
                    matrix22.matr[i - 1][j] = matrix_for_reverse.matr[i][j];
                }
            }
            else if (j > 2) {
                if (i < 2) {
                    matrix22.matr[i][j - 1] = matrix_for_reverse.matr[i][j];
                }
                else if (i > 2) {
                    matrix22.matr[i - 1][j - 1] = matrix_for_reverse.matr[i][j];
                }
            }
        }
    }
    //   matrix22.print_string();
    result.matr[2][2] = matrixDet(matrix22);

    Matrix<3, 3> matrix23;
    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 3; ++j) {
            if (i < 2) {
                matrix23.matr[i][j] = matrix_for_reverse.matr[i][j];
            }
            else if (i > 2) {
                matrix23.matr[i - 1][j] = matrix_for_reverse.matr[i][j];
            }
        }
    }
    // matrix23.print_string();
    result.matr[2][3] = -1 * matrixDet(matrix23);


    Matrix<3, 3> matrix30;
    for (int i = 0; i < 3; ++i) {
        for (int j = 1; j < 4; ++j) {
            matrix30.matr[i][j - 1] = matrix_for_reverse.matr[i][j];
        }
    }
    //  matrix30.print_string();
    result.matr[3][0] = -1 * matrixDet(matrix30);

    Matrix<3, 3> matrix31;
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 4; ++j) {
            if (j < 1) {
                matrix31.matr[i][j] = matrix_for_reverse.matr[i][j];
            }
            else if (j > 1) {
                matrix31.matr[i][j - 1] = matrix_for_reverse.matr[i][j];
            }
        }
    }
    // matrix31.print_string();
    result.matr[3][1] = matrixDet(matrix31);

    Matrix<3, 3> matrix32;
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 4; ++j) {
            if (j < 2) {
                matrix32.matr[i][j] = matrix_for_reverse.matr[i][j];
            }
            else if (j > 2) {
                matrix32.matr[i][j - 1] = matrix_for_reverse.matr[i][j];
            }
        }
    }
    // matrix32.print_string();
    result.matr[3][2] = -1 * matrixDet(matrix32);

    Matrix<3, 3> matrix33;
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            matrix33.matr[i][j] = matrix_for_reverse.matr[i][j];
        }
    }
    // matrix33.print_string();
    result.matr[3][3] = matrixDet(matrix33);

    //  result.print_string();

    cout << "transposed matrix " << endl;
    result = transposeMatrix(result);
    result.print_string();

    cout << "determ " << determinant << endl;

    if (determinant != 0) {
        result = multipleMatrix(1 / determinant, result);
    }
    else {
        result = multipleMatrix(1000, result);
    }


    return result;
}


template <int N>
Matrix<N, N> i_matrix() {
    Matrix<N, N> i_matr;
    for (int j = 0; j < N; ++j) {
        for (int k = 0; k < N; ++k) {
            if (j == k) {
                i_matr.matr[j][j] = 1;
            }
        }
    }

    return i_matr;
}

State<4> update_step_state(Matrix<4, 1> predicted_state, Matrix<2, 1> measured_state, Matrix<2, 1> prev_state,
    Matrix<4, 4> prior_predict_covar, Matrix<4, 4> measure_covar) {

    cout << "update step start " << endl;
    Matrix<4, 1> update_step_state_matrix;
    Matrix<4, 4> update_step_covar_matrix;

    cout << "measured state vector " << endl;
    Matrix<4, 1> measured_state_vector;     //Zk
    measured_state_vector.matr[0][0] = measured_state.matr[0][0];
    measured_state_vector.matr[1][0] = measured_state.matr[1][0];
    measured_state_vector.matr[2][0] = (measured_state.matr[0][0] - prev_state.matr[0][0]) / 0.1; ///actually 0.1 its timestamp and need to be a variable
    measured_state_vector.matr[3][0] = (measured_state.matr[1][0] - prev_state.matr[1][0]) / 0.1; // but for now we will leave it just 0.1
    measured_state_vector.print_string();

    cout << "innovation vector" << endl;
    //Vi(Xk)
    Matrix<4, 1> inovation_vector = minusMatrix(measured_state_vector, predicted_state);
    inovation_vector.print_string();

    Matrix<4, 4> s_matrix;

    s_matrix = multipleMatrix(multipleMatrix(getMeasurmentModelMatrix(), prior_predict_covar), transposeMatrix(getMeasurmentModelMatrix()));

    Matrix<4, 4> s_matrix_part_R = measure_covar;   //maybe todo

    s_matrix = addMatrix(s_matrix, s_matrix_part_R);
    cout << "finish S matrix" << endl;
    s_matrix.print_string();

    cout << "before K matrix" << endl;
    Matrix<4, 4> K_matrix = multipleMatrix(multipleMatrix(prior_predict_covar, transposeMatrix(getMeasurmentModelMatrix())),
        reverse_matrix(s_matrix));
    cout << "finish K matrix " << endl;
    K_matrix.print_string();

    update_step_state_matrix = addMatrix(measured_state_vector, multipleMatrix(K_matrix, inovation_vector));

    update_step_covar_matrix = addMatrix(i_matrix<4>(), multipleMatrix(K_matrix, getMeasurmentModelMatrix()));


    return { update_step_state_matrix, update_step_covar_matrix };

}


gausian func_update(long double prior_mean, long double prior_sigma2, long double measured_mean, long double measured_sigma2);

gausian func_Predict(long double prior_mean, long double prior_sigma2, long double dist, long double uncertenty2);

State<4> predict_step_state(Matrix<4, 1> start_state_vector, Matrix<4, 4> startPredictCovarMatrix, Matrix<4, 4> motion_matrix_01,
    Matrix<4, 4> Q_matrix_t_0_1) {
    cout << "prediction for next step" << endl;
    Matrix<4, 1> prediction_for_next_step = prediction_step_state_vector(start_state_vector, motion_matrix_01);
    prediction_for_next_step.print_string();

    cout << "covariance matrix predict" << endl;
    Matrix<4, 4> prediction_covar_matr = predict_step_covar_matrix(motion_matrix_01, startPredictCovarMatrix, Q_matrix_t_0_1);
    prediction_covar_matr.print_string();
    return { prediction_for_next_step, prediction_covar_matr };
}

void predict_and_update(State<4>& state_after_update, State<4>& predicted_state_next, Matrix<4, 4> motion_matrix_01,
    Matrix<4, 4> Q_matrix_t_0_1, int starti, int finishi, const double mes_r[500], const double mes_O[500]) {
    Matrix<2, 1> measured_state, prev_state;

    for (int i = starti; i < finishi; ++i) {
        cout << "prediction for step " << i << endl;

        predicted_state_next = predict_step_state(state_after_update.state_vector,
            state_after_update.covariance_matrix, motion_matrix_01,
            Q_matrix_t_0_1);

        SDL_SetRenderDrawColor(renderer, 255, 0, 0, SDL_ALPHA_OPAQUE);
        SDL_RenderDrawPoint(renderer, predicted_state_next.state_vector.matr[0][0] * 2 + 10, 210 - predicted_state_next.state_vector.matr[1][0] * 2);
        SDL_RenderPresent(renderer);
    	
        cout << "next prediction " << endl;
        predicted_state_next.state_vector.print_string();
        cout << "next prediction covar " << endl;
        predicted_state_next.covariance_matrix.print_string();

        measured_state.matr[0][0] = mes_r[i];
        measured_state.matr[1][0] = mes_O[i];


        this_thread::sleep_for(chrono::milliseconds(100));
        SDL_SetRenderDrawColor(renderer, 0, 100, 0, SDL_ALPHA_OPAQUE);
        SDL_RenderDrawPoint(renderer, measured_state.matr[0][0] * 2 + 10, 210 - measured_state.matr[1][0] * 2);
        SDL_RenderPresent(renderer);

        prev_state.matr[0][0] = state_after_update.state_vector.matr[0][0];
        prev_state.matr[1][0] = state_after_update.state_vector.matr[1][0];

        cout << "update for step " << i << endl;

        state_after_update = update_step_state(predicted_state_next.state_vector, measured_state, prev_state,
            predicted_state_next.covariance_matrix, motion_matrix_01);

        cout << "updated vector " << endl;
        state_after_update.state_vector.print_string();
        cout << "updated covar " << endl;
        state_after_update.covariance_matrix.print_string();

        cout << "go to next step " << endl;
    }
}

void first_step_state_vector_calculation(Matrix<4, 1>& start_state_vector, double x, double y) {
    start_state_vector.matr[0][0] = x;
    start_state_vector.matr[1][0] = y;
    start_state_vector.matr[2][0] = 0;
    start_state_vector.matr[3][0] = 0;
}

void first_step_state_vector_calculation(Matrix<4, 1>& start_state_vector, double x, double y, double nextx, double nexty, double dt) {
    start_state_vector.matr[0][0] = nextx;
    start_state_vector.matr[1][0] = nexty;
    start_state_vector.matr[2][0] = (nextx - x) / dt;
    start_state_vector.matr[3][0] = (nexty - y) / dt;
}

//if we will want start from no measured case, just pass 0, 0 but for EKF its really dangerous
void first_state_using_1_measurment(Matrix<4, 1>& start_state_vector, Matrix<4, 4> startPredictCovarMatrix, double x, double y) {

    first_step_state_vector_calculation(start_state_vector, x, y);

    SDL_SetRenderDrawColor(renderer, 0, 100, 0, SDL_ALPHA_OPAQUE);
    SDL_RenderDrawPoint(renderer, x*2 + 10, 210 - y*2);
    SDL_RenderPresent(renderer);

    startPredictCovarMatrix.matr[0][0] = R_STANDART_DEVIATION * (1 + cos(PSI_STANDART_DEVIATION)); //equal mesurment error since we use first megurment for start value
    startPredictCovarMatrix.matr[1][1] = R_STANDART_DEVIATION * (1 + sin(PSI_STANDART_DEVIATION)); //i will use just diagonal values for start aproximation
    startPredictCovarMatrix.matr[2][2] = start_state_vector.matr[0][0];
    startPredictCovarMatrix.matr[3][3] = start_state_vector.matr[1][0];
    cout << "start prediction covar matrix " << endl;
    startPredictCovarMatrix.print_string();

    //   return {start_state_vector, startPredictCovarMatrix};
}

void first_state_using_2(Matrix<4, 1>& start_state_vector, Matrix<4, 4> startPredictCovarMatrix, double x, double y, double nextx, double nexty) {

    first_step_state_vector_calculation(start_state_vector, x, y, nextx, nexty, 0.1);

    startPredictCovarMatrix.matr[0][0] = R_STANDART_DEVIATION * (1 + cos(PSI_STANDART_DEVIATION)); //equal mesurment error since we use first megurment for start value
    startPredictCovarMatrix.matr[1][1] = R_STANDART_DEVIATION * (1 + sin(PSI_STANDART_DEVIATION)); //i will use just diagonal values for start aproximation
    startPredictCovarMatrix.matr[2][2] = start_state_vector.matr[0][0];
    startPredictCovarMatrix.matr[3][3] = start_state_vector.matr[1][0];
    cout << "start prediction covar matrix " << endl;
    startPredictCovarMatrix.print_string();

    //   return {start_state_vector, startPredictCovarMatrix};
}

void data_for_first_update(Matrix<2, 1>& measured_state, Matrix<2, 1>& prev_state, double first_x, double first_y, double next_x, double next_y) {
    measured_state.matr[0][0] = first_x;
    measured_state.matr[1][0] = first_y;

    ///probably for start it is
    prev_state.matr[0][0] = next_x;
    prev_state.matr[1][0] = next_y;

}

//for motion data
gausian func_Predict(long double prior_mean, long double prior_sigma2, long double dist, long double uncertenty2) {
    long double new_mean = prior_mean + dist;
    long double new_variance = prior_sigma2 + uncertenty2;
    return { new_mean, new_variance };
}

void renderText(TTF_Font* font, const char* textStr, SDL_Color color, int x, int y) {
    SDL_Surface* text;
    text = TTF_RenderText_Solid(font, textStr, color);
    if (!text)
    {
        cout << "Failed to render text: " << TTF_GetError() << endl;
    }
    SDL_Texture* text_texture;
    text_texture = SDL_CreateTextureFromSurface(renderer, text);
    SDL_Rect dest = { x, y, text->w, text->h };
    SDL_RenderCopy(renderer, text_texture, NULL, &dest);
}

int main()
{
	SDL_SetMainReady();
	cout << "Hello!" << endl;
	

	const double INIT_SIGMA_R_2 = pow(R_STANDART_DEVIATION, 2);
	const double INIT_SIGMA_PSI_2 = pow(PSI_STANDART_DEVIATION, 2);

	double mes_r[500] = {
		50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77,
		78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 103, 104,
		105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115, 116, 117, 118, 119, 120, 121, 122, 123, 124, 125, 126,
		127, 128, 129, 130, 131, 132, 133, 134, 135, 136, 137, 138, 139, 140, 141, 142, 143, 144, 145, 146, 147, 148,
		149, 149.999987307656, 150.999898462407, 151.999657317147, 152.999187736368, 153.998413601960, 154.997258819004,
		155.995647321570, 156.993503078507, 157.990750099237, 158.987312439538, 159.983114207328, 160.978079568450,
		161.972132752437, 162.965198058291, 163.957199860245, 164.948062613520, 165.937710860081, 166.926069234384,
		167.913062469111, 168.898615400906, 169.882652976098, 170.865100256415, 171.845882424693, 172.824924790573,
		173.802152796188, 174.777492021843, 175.750868191679, 176.722207179331, 177.691435013576, 178.658477883961,
		179.623262146428, 180.585714328920, 181.545761136976, 182.503329459316, 183.458346373404, 184.410739151005,
		185.360435263721, 186.307362388517, 187.251448413225, 188.192621442041, 189.130809800993, 190.065942043406,
		190.997946955339, 191.926753561010, 192.852291128200, 193.774489173639, 194.693277468376, 195.608586043125,
		196.520345193594, 197.428485485795, 198.332937761327, 199.233633142648, 200.130503038318, 201.023479148222,
		201.912493468772, 202.797478298085, 203.678366241141, 204.555090214911, 205.427583453471, 206.295779513082,
		207.163975572694, 208.036468811254, 208.913192785024, 209.794080728080, 210.679065557393, 211.568079877942,
		212.461055987846, 213.357925883517, 214.258621264838, 215.163073540370, 216.071213832571, 216.982972983040,
		217.898281557789, 218.817069852526, 219.739267897965, 220.664805465155, 221.593612070826, 222.525616982759,
		223.460749225172, 224.398937584124, 225.340110612939, 226.284196637648, 227.231123762444, 228.180819875160,
		229.133212652761, 230.088229566849, 231.045797889189, 232.005844697245, 232.968296879737, 233.933081142204,
		234.900124012589, 235.869351846834, 236.840690834486, 237.814067004322, 238.789406229976, 239.766634235592,
		240.745676601472, 241.726458769750, 242.708906050067, 243.692943625259, 244.678496557054, 245.665489791781,
		246.653848166083, 247.643496412645, 248.634359165920, 249.626360967874, 250.619426273728, 251.613479457715,
		252.608444818836, 253.604246586627, 254.600808926928, 255.598055947658, 256.595911704595, 257.594300207161,
		258.593145424205, 259.592371289797, 260.591901709018, 261.591660563758, 262.591571718509, 263.591559026165,
		264.591559026165, 265.591559026165, 266.591559026165, 267.591559026165, 268.591559026165, 269.591559026165,
		270.591559026165, 271.591559026165, 272.591559026165, 273.591559026165, 274.591559026165, 275.591559026165,
		276.591559026165, 277.591559026165, 278.591559026165, 279.591559026165, 280.591559026165, 281.591559026165,
		282.591559026165, 283.591559026165, 284.591559026165, 285.591559026165, 286.591559026165, 287.591559026165,
		288.591559026165, 289.591559026165, 290.591559026165, 291.591559026165, 292.591559026165, 293.591559026165,
		294.591559026165, 295.591559026165, 296.591559026165, 297.591559026165, 298.591559026165, 299.591559026165,
		300.591559026165, 301.591559026165, 302.591559026165, 303.591559026165, 304.591559026165, 305.591559026165,
		306.591559026165, 307.591559026165, 308.591559026165, 309.591559026165, 310.591559026165, 311.591559026165,
		312.591559026165, 313.591559026165, 314.591559026165, 315.591559026165, 316.591559026165, 317.591559026165,
		318.591559026165, 319.591559026165, 320.591559026165, 321.591559026165, 322.591559026165, 323.591559026165,
		324.591559026165, 325.591559026165, 326.591559026165, 327.591559026165, 328.591559026165, 329.591559026165,
		330.591559026165, 331.591559026165, 332.591559026165, 333.591559026165, 334.591559026165, 335.591559026165,
		336.591559026165, 337.591559026165, 338.591559026165, 339.591559026165, 340.591559026165, 341.591559026165,
		342.591559026165, 343.591559026165, 344.591554456910, 345.591522472275, 346.591435657632, 347.591266599855,
		348.590987888076, 349.590572114428, 350.589991874806, 351.589219769613, 352.588228404510, 353.586990391172,
		354.585478348033, 355.583664901043, 356.581522684412, 357.579024341366, 358.576142524893, 359.572849898493,
		360.569119136932, 361.564922926984, 362.560233968185, 363.555024973581, 364.549268670473, 365.542937801169,
		366.536005123727, 367.528443412705, 368.520225459906, 369.511324075124, 370.501712086889, 371.491362343216,
		372.480247712340, 373.468341083472, 374.455615367531, 375.442043497894, 376.427598431137, 377.412253147772,
		378.395980652994, 379.378753977415, 380.360546177809, 381.341330337847, 382.321079568835, 383.299767010453,
		384.277365831491, 385.253849230583, 386.229190436943, 387.203362711100, 388.176339345629, 389.148093665883,
		390.118599030728, 391.087828833268, 392.055756501578, 393.022355499433, 393.987599327033, 394.951461521731,
		395.913915658758, 396.874935351947, 397.834494254458, 398.792566059500, 399.749124501049, 400.704143354572,
		401.657596437746, 402.609457611173, 403.559700779097, 404.508299890122, 405.455228937925, 406.400461961967,
		407.343973048207, 408.285736329812, 409.225725987867, 410.163916252080, 411.100281401490, 412.034795765175,
		412.967433722950, 413.898169706076, 414.826978197955, 415.753833734833, 416.678710906497, 417.601584356973,
		418.522428785219, 419.441218945819, 420.357929649677, 421.272535764706, 422.185012216518, 423.095333989109,
		424.003476125548, 424.909413728661, 425.813121961709, 426.714576049076, 427.613751276941, 428.510622993963,
		429.405166611949, 430.297357606535, 431.187171517854, 432.074583951208, 432.959570577736, 433.842107135084,
		434.722169428065, 435.599733329329, 436.474774780016, 437.347269790423, 438.217194440660, 439.084524881302,
		439.951855321944, 440.821779972180, 441.694274982588, 442.569316433275, 443.446880334538, 444.326942627520,
		445.209479184868, 446.094465811396, 446.981878244750, 447.871692156069, 448.763883150655, 449.658426768641,
		450.555298485662, 451.454473713528, 452.355927800895, 453.259636033943, 454.165573637056, 455.073715773495,
		455.984037546086, 456.896513997898, 457.811120112927, 458.727830816785, 459.646620977385, 460.567465405631,
		461.490338856106, 462.415216027771, 463.342071564649, 464.270880056528, 465.201616039654, 466.134253997429,
		467.068768361114, 468.005133510524, 468.943323774737, 469.883313432791, 470.825076714397, 471.768587800637,
		472.713820824679, 473.660749872481, 474.609348983507, 475.559592151431, 476.511453324858, 477.464906408032,
		478.419925261555, 479.376483703104, 480.334555508146, 481.294114410657, 482.255134103846, 483.217588240873,
		484.181450435570, 485.146694263171, 486.113293261026, 487.081220929336, 488.050450731876, 489.020956096721,
		489.992710416975, 490.965687051504, 491.939859325661, 492.915200532021, 493.891683931113, 494.869282752151,
		495.847970193769, 496.827719424757, 497.808503584794, 498.790295785188, 499.773069109610, 500.756796614831,
		501.741451331467, 502.727006264710, 503.713434395073, 504.700708679132, 505.688802050264, 506.677687419388,
		507.667337675714, 508.657725687480, 509.648824302698, 510.640606349899, 511.633044638877, 512.626111961435,
		513.619781092131, 514.614024789023, 515.608815794419, 516.604126835620, 517.599930625672, 518.596199864110,
		519.592907237711, 520.590025421238, 521.587527078191, 522.585384861561, 523.583571414571, 524.582059371432,
		525.580821358094, 526.579829992991, 527.579057887797, 528.578477648176, 529.578061874528, 530.577783162748,
		531.577614104972, 532.577527290329, 533.577495305694, 534.577490736439
	};

	typedef numeric_limits<double> dbl;

	cout.precision(13);
	cout << "last r = " << mes_r[499] << endl;

	double mes_O[500] = {
		50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50,
		50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50,
		50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50,
		50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50.0043632954396, 50.0174528494764,
		50.0392676652902, 50.0698060815984, 50.1090657727824, 50.1570437490646, 50.2137363567363, 50.2791392784361,
		50.3532475334782, 50.4360554782321, 50.5275568065522, 50.6277445502580, 50.7366110796650, 50.8541481041655,
		50.9803466728600, 51.1151971752389, 51.2586893419144, 51.4108122454023, 51.5715543009547, 51.7409032674416,
		51.9188462482837, 52.1053696924339, 52.3004593954099, 52.5041005003754, 52.7162774992719, 52.9369742339995,
		53.1661738976473, 53.4038590357735, 53.6500115477347, 53.9046126880640, 54.1676430678989, 54.4390826564577,
		54.7189107825646, 55.0071061362245, 55.3036467702454, 55.6085101019100, 55.9216729146953, 56.2431113600404,
		56.5728009591634, 56.9107166049249, 57.2568325637400, 57.6111224775386, 57.9735593657719, 58.3441156274678,
		58.7227630433321, 59.1094727778983, 59.5042153817228, 59.9069607936279, 60.3176783429914, 60.7363367520817,
		61.1629041384400, 61.5973480173083, 62.0396353041032, 62.4897323169355, 62.9476047791750, 63.4132178220609,
		63.8865359873574, 64.3675232300537, 64.8561429211090, 65.3523578502421, 65.8485727793751, 66.3371924704304,
		66.8181797131267, 67.2914978784232, 67.7571109213092, 68.2149833835487, 68.6650803963809, 69.1073676831759,
		69.5418115620442, 69.9683789484024, 70.3870373574927, 70.7977549068562, 71.2005003187614, 71.5952429225859,
		71.9819526571520, 72.3606000730164, 72.7311563347122, 73.0935932229456, 73.4478831367441, 73.7939990955593,
		74.1319147413207, 74.4616043404437, 74.7830427857889, 75.0962055985741, 75.4010689302387, 75.6976095642597,
		75.9858049179196, 76.2656330440265, 76.5370726325852, 76.8001030124201, 77.0547041527495, 77.3008566647106,
		77.5385418028369, 77.7677414664847, 77.9884382012122, 78.2006152001087, 78.4042563050743, 78.5993460080502,
		78.7858694522005, 78.9638124330425, 79.1331613995295, 79.2939034550818, 79.4460263585698, 79.5895185252452,
		79.7243690276241, 79.8505675963186, 79.9681046208191, 80.0769711502261, 80.1771588939320, 80.2686602222520,
		80.3514681670059, 80.4255764220480, 80.4909793437478, 80.5476719514196, 80.5956499277018, 80.6349096188857,
		80.6654480351940, 80.6872628510078, 80.7003524050446, 80.7047157004842, 80.7047157004842, 80.7047157004842,
		80.7047157004842, 80.7047157004842, 80.7047157004842, 80.7047157004842, 80.7047157004842, 80.7047157004842,
		80.7047157004842, 80.7047157004842, 80.7047157004842, 80.7047157004842, 80.7047157004842, 80.7047157004842,
		80.7047157004842, 80.7047157004842, 80.7047157004842, 80.7047157004842, 80.7047157004842, 80.7047157004842,
		80.7047157004842, 80.7047157004842, 80.7047157004842, 80.7047157004842, 80.7047157004842, 80.7047157004842,
		80.7047157004842, 80.7047157004842, 80.7047157004842, 80.7047157004842, 80.7047157004842, 80.7047157004842,
		80.7047157004842, 80.7047157004842, 80.7047157004842, 80.7047157004842, 80.7047157004842, 80.7047157004842,
		80.7047157004842, 80.7047157004842, 80.7047157004842, 80.7047157004842, 80.7047157004842, 80.7047157004842,
		80.7047157004842, 80.7047157004842, 80.7047157004842, 80.7047157004842, 80.7047157004842, 80.7047157004842,
		80.7047157004842, 80.7047157004842, 80.7047157004842, 80.7047157004842, 80.7047157004842, 80.7047157004842,
		80.7047157004842, 80.7047157004842, 80.7047157004842, 80.7047157004842, 80.7047157004842, 80.7047157004842,
		80.7047157004842, 80.7047157004842, 80.7047157004842, 80.7047157004842, 80.7047157004842, 80.7047157004842,
		80.7047157004842, 80.7047157004842, 80.7047157004842, 80.7047157004842, 80.7047157004842, 80.7047157004842,
		80.7047157004842, 80.7047157004842, 80.7047157004842, 80.7047157004842, 80.7047157004842, 80.7047157004842,
		80.7020977125873, 80.6942438206702, 80.6811542400513, 80.6628293295881, 80.6392695916671, 80.6104756721906,
		80.5764483605584, 80.5371885896464, 80.4926974357811, 80.4429761187099, 80.3880260015679, 80.3278485908402,
		80.2624455363210, 80.1918186310681, 80.1159698113538, 80.0349011566117, 79.9486148893800, 79.8571133752402,
		79.7603991227527, 79.6584747833875, 79.5513431514518, 79.4390071640135, 79.3214699008204, 79.1987345842159,
		79.0708045790508, 78.9376833925907, 78.7993746744203, 78.6558822163429, 78.5072099522767, 78.3533619581470,
		78.1943424517742, 78.0301557927584, 77.8608064823598, 77.6862991633754, 77.5066386200113, 77.3218297777522,
		77.1318777032259, 76.9367876040644, 76.7365648287612, 76.5312148665251, 76.3207433471288, 76.1051560407555,
		75.8844588578399, 75.6586578489070, 75.4277592044053, 75.1917692545380, 74.9506944690886, 74.7045414572442,
		74.4533169674139, 74.1970278870441, 73.9356812424293, 73.6692841985199, 73.3978440587253, 73.1213682647141,
		72.8398643962098, 72.5533401707833, 72.2618034436409, 71.9652622074092, 71.6637245919160, 71.3571988639674,
		71.0456934271210, 70.7292168214558, 70.4077777233377, 70.0813849451821, 69.7500474352117, 69.4137742772119,
		69.0725746902810, 68.7264580285782, 68.3754337810666, 68.0195115712531, 67.6587011569251, 67.2930124298822,
		66.9224554156657, 66.5470402732832, 66.1667772949306, 65.7816769057094, 65.3917496633414, 64.9970062578790,
		64.5974575114119, 64.1931143777709, 63.7839879422272, 63.3700894211887, 62.9514301618926, 62.5280216420938,
		62.0998754697508, 61.6670033827072, 61.2294172483699, 60.7871290633839, 60.3401509533031, 59.8884951722583,
		59.4321741026208, 58.9712002546632, 58.5055862662165, 58.0353449023234, 57.5604890548883, 57.0810317423243,
		56.5969861091957, 56.1083654258581, 55.6151830880943, 55.1174526167473, 54.6197221454003, 54.1265398076365,
		53.6379191242989, 53.1538734911703, 52.6744161786063, 52.1995603311712, 51.7293189672781, 51.2637049788314,
		50.8027311308738, 50.3464100612363, 49.8947542801915, 49.4477761701107, 49.0054879851247, 48.5679018507874,
		48.1350297637438, 47.7068835914008, 47.2834750716020, 46.8648158123059, 46.4509172912674, 46.0417908557238,
		45.6374477220827, 45.2378989756156, 44.8431555701532, 44.4532283277852, 44.0681279385641, 43.6878649602114,
		43.3124498178289, 42.9418928036124, 42.5762040765695, 42.2153936622415, 41.8594714524280, 41.5084472049164,
		41.1623305432136, 40.8211309562827, 40.4848577982829, 40.1535202883125, 39.8271275101569, 39.5056884120388,
		39.1892118063736, 38.8777063695272, 38.5711806415786, 38.2696430260854, 37.9731017898537, 37.6815650627113,
		37.3950408372847, 37.1135369687805, 36.8370611747693, 36.5656210349747, 36.2992239910653, 36.0378773464505,
		35.7815882660807, 35.5303637762504, 35.2842107644060, 35.0431359789566, 34.8071460290892, 34.5762473845876,
		34.3504463756547, 34.1297491927391, 33.9141618863658, 33.7036903669695, 33.4983404047333, 33.2981176294302,
		33.1030275302687, 32.9130754557423, 32.7282666134833, 32.5486060701192, 32.3740987511348, 32.2047494407362,
		32.0405627817204, 31.8815432753476, 31.7276952812179, 31.5790230171517, 31.4355305590743, 31.2972218409039,
		31.1641006544438, 31.0361706492787, 30.9134353326742, 30.7958980694811, 30.6835620820428, 30.5764304501071,
		30.4745061107419, 30.3777918582544, 30.2862903441147, 30.2000040768829, 30.1189354221408, 30.0430866024265,
		29.9724596971736, 29.9070566426544, 29.8468792319268, 29.7919291147847, 29.7422077977135, 29.6977166438482,
		29.6584568729362, 29.6244295613040, 29.5956356418275, 29.5720759039066, 29.5537509934433, 29.5406614128244,
		29.5328075209073, 29.5301895330105
	};

	cout << "last second " << mes_O[499] << endl;
	/*

	    long double mu, x;
	    long double sigma;

	    long double fx = func_fx(10.0, 2.0, 8.0);

	    cout << "double fx = " << fx << endl;
	*/

	/// <summary>
	/// hear I will start to invoke my window with plot
	/// </summary>
	/// <returns></returns>


	double* max_r = max_element(mes_r, mes_r + 500);
	double* max_O = max_element(mes_O, mes_O + 500);
	cout << "max x" << *max_r << " max y " << *max_O << endl;
	double max_X = *max_r;
	double max_Y = *max_O; //TODO use than programmable for render margin but for now it just hard code//

	SDL_Window* window; // Declare a pointer

	SDL_Init(SDL_INIT_VIDEO); // Initialize SDL2
    if (TTF_Init() < 0)
    {
        cout << "Error initializing SDL_ttf: " << TTF_GetError() << endl;
    }
	TTF_Font* font;

	font = TTF_OpenFont("D:\\roboto\\Roboto-Black.ttf", 14);
	if (!font)
	{
		cout << "Failed to load font: " << TTF_GetError() << endl;
	}

	// Create an application window with the following settings:
	window = SDL_CreateWindow(
		"EKF real-time plot", // window title
		SDL_WINDOWPOS_UNDEFINED, // initial x position
		SDL_WINDOWPOS_UNDEFINED, // initial y position
		1280, // width, in pixels
		720, // height, in pixels
		SDL_WINDOW_OPENGL // flags - see below
	);

	// Check that the window was successfully created
	if (window == NULL)
	{
		// In the case that the window could not be made...
		printf("Could not create window: %s\n", SDL_GetError());
	}

	renderer = SDL_CreateRenderer(window, 0, SDL_WINDOW_OPENGL);


	// The window is open: could enter program loop here (see SDL_PollEvent())
	SDL_Event event;

	SDL_SetRenderDrawColor(renderer, 255, 255, 255, SDL_ALPHA_OPAQUE);
	SDL_RenderClear(renderer);

	SDL_SetRenderDrawColor(renderer, 0, 0, 0, SDL_ALPHA_OPAQUE);
	SDL_RenderDrawLine(renderer, 10, 10, 10, 210);
	SDL_RenderDrawLine(renderer, 10, 10, 1110, 10);
	SDL_RenderDrawLine(renderer, 1110, 10, 1110, 210);
	SDL_RenderDrawLine(renderer, 10, 210, 1110, 210);

    SDL_RenderDrawLine(renderer, 10, 110, 20, 110);
	SDL_RenderDrawLine(renderer, 10, 60, 20, 60);
	SDL_RenderDrawLine(renderer, 10, 160, 20, 160);

	SDL_RenderDrawLine(renderer, 1010, 200, 1010, 210);
	SDL_RenderDrawLine(renderer, 510, 200, 510, 210);
	SDL_RenderDrawLine(renderer, 260, 200, 260, 210);
	SDL_RenderDrawLine(renderer, 760, 200, 760, 210);

	SDL_RenderPresent(renderer);
	//  SDL_Delay(1000);

    // Set color to black
    SDL_Color color = { 0, 0, 0 };
   
    renderText(font, "0", color, 10, 210);
    renderText(font, "125", color, 260, 210);
    renderText(font, "250", color, 510, 210);
    renderText(font, "375", color, 760, 210);
    renderText(font, "500", color, 1010, 210);

    renderText(font, "25", color, 12, 160);
    renderText(font, "50", color, 12, 110);
    renderText(font, "75", color, 12, 60);


	Matrix<4, 4> meagurment_covar_temp; //TODO
	meagurment_covar_temp.matr[0][0] = R_STANDART_DEVIATION * (1 + cos(PSI_STANDART_DEVIATION));
	meagurment_covar_temp.matr[1][1] = R_STANDART_DEVIATION * (1 + sin(PSI_STANDART_DEVIATION));
	cout << "measurment covar matrix " << endl;
	meagurment_covar_temp.print_string();

	Matrix<4, 4> motion_matrix_01 = get_motion_matrix(0.1);

	cout << "motion covar matrix" << endl;
	motion_matrix_01.print_string();
	/*  getVelocityNoiseCovarianceMatrix(1,0.1).print_string();

	  get_diff_t_velocity_covariance_matrix(1, 0.1).print_string();
  */

	Matrix<4, 4> Q_matrix_t_0_1 = get_velocity_noise_covariance_matrix(0.3, 0.1);
	//TODO q lowercase that we will change for some p. of task

	cout << "Q matrix " << endl;
	Q_matrix_t_0_1.print_string();

	Matrix<4, 1> start_state_vector;
	Matrix<4, 4> startPredictCovarMatrix;
	Matrix<2, 1> measured_state;
	Matrix<2, 1> prev_state;
	State<4> predicted_state;
	State<4> state_after_update;

	///first step
	first_state_using_1_measurment(start_state_vector, startPredictCovarMatrix, mes_r[0], mes_O[0]);

	predicted_state = predict_step_state(start_state_vector, startPredictCovarMatrix, motion_matrix_01, Q_matrix_t_0_1);


	SDL_SetRenderDrawColor(renderer, 255, 0, 0, SDL_ALPHA_OPAQUE);
	SDL_RenderDrawPoint(renderer, predicted_state.state_vector.matr[0][0] * 2 + 10,
	                    210 - predicted_state.state_vector.matr[1][0] * 2);
	SDL_RenderPresent(renderer);

	data_for_first_update(measured_state, prev_state, mes_r[0], mes_O[0], mes_r[1], mes_O[1]);

	SDL_SetRenderDrawColor(renderer, 0, 100, 0, SDL_ALPHA_OPAQUE);
	SDL_RenderDrawPoint(renderer, mes_r[1] * 2 + 10, 210 - mes_O[1] * 2);
	SDL_RenderPresent(renderer);

	state_after_update = update_step_state(predicted_state.state_vector, measured_state, prev_state,
	                                       predicted_state.covariance_matrix, meagurment_covar_temp);

	cout << "state after update" << endl;
	state_after_update.state_vector.print_string();
	cout << "covar after predict " << endl;
	state_after_update.covariance_matrix.print_string();

	///for the all next steps

	predict_and_update(state_after_update, predicted_state, motion_matrix_01, Q_matrix_t_0_1, 2, 499, mes_r, mes_O);


	cin.get();
	//Close winwow with graphics on enter
	// Close and destroy the window
  //  SDL_DestroyTexture(text_texture0);
  //  SDL_FreeSurface(text0);
	if (renderer)
	{
		SDL_DestroyRenderer(renderer);
	}
	if (window)
	{
		SDL_DestroyWindow(window);
	}
	// Clean up
	SDL_Quit();

	return 0;
}

long double func_new_sigma2(long double prior_sigma2, long double measured_sigma2) {
    return 1 / (1 / prior_sigma2 + 1 / measured_sigma2);
}

long double func_new_mean(long double prior_mean, long double prior_sigma2, long double measured_mean, long double measured_sigma2) {
    return (prior_mean * measured_sigma2 + prior_sigma2 * measured_mean) / (prior_sigma2 + measured_sigma2);
}

long double func_fx(long double mu, long double sigma2, long double x) {


    // long double sigma2 = pow(sigma, 2);
    long double result = (1 / sqrt(2 * 3.1415926535 * sigma2)) * exp(-0.5 * pow((x - mu), 2) / sigma2);

    return result;
}


gausian func_update(long double prior_mean, long double prior_sigma2, long double measured_mean, long double measured_sigma2) {

    //after update
    long double posterior_mean = func_new_mean(prior_mean, prior_sigma2, measured_mean, measured_sigma2);   //new mean
    long double posterior_sigma2 = func_new_sigma2(prior_sigma2, measured_sigma2);  //mew variance

    gausian gaus(posterior_mean, posterior_sigma2);
    return gaus;

}