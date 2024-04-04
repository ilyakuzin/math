import java.util.Arrays;

public class Main {
    private static final int M = 5; // разбиение по оси t
    private static final int N = 100; // разбиение по оси x.
    private static final double HALF_ONE = 0.5;
    private static final double T = 1.0;
    private static final double h = 1.0 / (N - 1);
    private static final double tau = T / (M - 1);
    private static final double hSquare = h * h;

    private static final double A = tau / hSquare;
    private static final double B = (2 * tau / hSquare + 3 / 2.);
    private static final double C = A;

    private static double mu_1(double t) {
        return Math.sin(t);
    }

    private static double mu_2(double t) {
        return Math.cos(t) * Math.sin(t);
    }

    private static double phi_1(double t) {
        return 0;
    }

    private static double phi_2(double x, double t) {
        return Math.cos(x) * Math.sin(t);
    }

    private static double f(double x, double t) {
        return Math.cos(x) * Math.sin(t) + Math.cos(x) * Math.cos(t);
    }

    private static double resultFunction(double x, double t) {
        return Math.cos(x) * Math.sin(t);
    } //правильное решение

    public static void main(String[] args) {

        double[] result_i = new double[N]; // значения функции-ответа на последнем слое (ответ).

        // Вычисляем y_i на самом нижнем слое. Заодно вычисляем результат на самом верхнем слое.
        double[] prevY_i = new double[N]; // y_(i-1).
        double[] y_i = new double[N];     // y_i
        double[] nextY_i = new double[N]; // y_(i+1).
        double x = 0.;


        //результирующая функция + у на 0 слое + краевые условия
        for (int i = 0; i < N; i++) { // j = 0.
            result_i[i] = resultFunction(x, T);
            prevY_i[i] = phi_1(x); //вычисление на 0 слое
            //y_i[i] = phi_1(x); //вычисление на 0 слое
            x += h;
        }
        prevY_i[0] = mu_1(tau); //краевые условия -- слева
        prevY_i[N - 1] = mu_2(tau); //справа


        double t = tau; // t^j.
        double[] D_i = new double[N];

        {
            double[] f_i = new double[N];     // f^0_i.
            double[] f_i_jpp = new double[N]; // f^(j+1)_i.
            double[] phi_2_i = new double[N]; //краевые условия
            x = 0; // x_i.

            for (int i = 1; i < N; i++) { // Прямой ход метода прогонки: вычисление f_i.
                f_i[i] = f(x, 0);
                f_i_jpp[i] = f(x, t);
                phi_2_i[i] = phi_2(x, t);
                x += h;
            }

            //y_i = new double[N];

            //вычисление U_1

            //y_i[N - 1] = mu_2(t + tau);
            x = 0;
            for (int i = 0; i < N - 1; i++) {
                y_i[i] = prevY_i[i] + tau * (phi_2_i[i] + f_i[i]); //вычислили значения на первом слое
                result_i[i] = resultFunction(x, T);

                x += h;
            }


            for (int i = 0; i < N; i++) { // Обновление коэффициентов D_i по y_i и y_(i-1).
                D_i[i] = (2 * y_i[i] - HALF_ONE * prevY_i[i] + tau * f_i_jpp[i]);
            }
        }

        for (int j = 2; j <= M; j++) { // Проход снизу вверх от слоя к слою.
            double[] f_i_jpp = new double[N]; // f^(j+1)_i.
            double[] alpha_i = new double[N];
            double[] beta_i = new double[N];
            alpha_i[0] = 0;
            beta_i[0] = mu_1(t);
            x = 0; // x_i.

            for (int i = 1; i < N - 1; i++) { // Прямой ход метода прогонки: вычисление f_i_jpp, alpha_i, beta_i.
                f_i_jpp[i] = f(x, t + tau);
                double denom = B - C * alpha_i[i - 1];
                System.out.println(alpha_i[i]);
                System.out.println();
                alpha_i[i] = A / denom; //по формуле
                beta_i[i] = (D_i[i - 1] + C * beta_i[i - 1]) / denom; //по формуле стр.19
                x += h;
            }

            nextY_i[N - 1] = mu_2(t);
            for (int i = N - 2; i >= 0; i--) { // Обратный ход метода прогонки: вычисление y_i по y_(i+1), alpha_(i+1) и beta_(i+1).
                nextY_i[i] = alpha_i[i] * nextY_i[i + 1] + beta_i[i];
            }
            x = 0.;

            //копируем массивы
            prevY_i = Arrays.copyOf(y_i, y_i.length);
            y_i = Arrays.copyOf(nextY_i, nextY_i.length);

            for (int i = 0; i < N; i++) { // Обновление коэффициентов D_i по y_i и y_(i-1).
                D_i[i] = (2 * y_i[i] - HALF_ONE * prevY_i[i] + tau * f_i_jpp[i]);
            }
            //обновляем значение t для вычислений
            t += tau;
        }


        double z = 0.;
        for (int i = 0; i < N; i++) {
            z = Math.max(z, Math.abs(y_i[i] - result_i[i]));
        }

        double result_norma = 0.;
        for (int i = 0; i < N; i++) {
            result_norma = Math.max(result_norma, Math.abs(result_i[i]));
            System.out.printf("|%f - %f| = %f \n", nextY_i[i], result_i[i], Math.abs(nextY_i[i] - result_i[i]));
        }

        System.out.println("Кубическая норма разности полученного решения и правильного, z: ");
        System.out.println(z);

        System.out.println("Отношение z к норме правильного решения: ");
        System.out.println(z / result_norma);

    }


}
