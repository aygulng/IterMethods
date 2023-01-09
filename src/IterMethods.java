public class IterMethods {
    static int n = 10;
    static double h = 1.0 / n;
    static double eps = Math.pow(10, -2);
    static double[] x = new double[n + 1];
    static double[] ai = new double[n + 1];
    static double[] qi = new double[n + 1];
    static double[] fi = new double[n + 1];
    static double[] alpha = new double[n + 1];
    static double[] beta = new double[n + 1];
    static double[] a_big = new double[n + 1];
    static double[] b_big = new double[n + 1];
    static double[] c_big = new double[n + 1];
    static double[][] A = new double[3] [n + 1];
    static double[] y = new double[n + 1];
    static double[] y_jacobi = new double[n + 1];
    static double[] y_zeid = new double[n + 1];
    static double[] y_relax = new double[n + 1];
    static double[] y_sin = new double[n + 1];
    static int k;
    static int best_k = 99999;
    static double best_omega;
    public static double Px(double x)
    {
        return 1 + x;
    }
    public static double Sin(double x)
    {
        return Math.sin(x);
    }
    public static double Qx(double x)
    {
        return x + 1;
    }
    public static double Ux(double x)
    {
        return x *  Math.pow((1 - x), 3);
    }
    public static double Fx(double x)
    {
        return ( - Math.pow(x, 5) + 2 * Math.pow(x, 4)
                + 16 * Math.pow(x, 3) - 17 * x * x - 5 * x + 5);
    }
    public static void outputMat(double[][] matrix) {
        for (double[] row : matrix) {
            for (double el : row)
                System.out.print(" " + el);
            System.out.println();
        }
    }
    public static void main(String[] args) {
        alpha[1] = 0; // так как c0=0
        beta[1] = 0; //так как f0=0
        y[0] = 0;
        y[n] = 0;
        // вычисление a[i],q[i],f[i]
        for (int i = 0; i <= n; i++)
        {
            ai[i] = Px(i * h);
            qi[i] = Qx(i * h);
            fi[i] = Fx(i * h) * Math.pow(h, 2);
        }
        for (int i = 0; i < n; i++)
            x[i] = i * h;

        // вычисление A,B,C
        for (int i = 1; i < n; i++)
        {
            a_big[i] = -ai[i];
            b_big[i] = -(ai[i] + ai[i + 1] + h * h * qi[i]);
            c_big[i] = -ai[i + 1];
        }
        for (int i = 0; i < n; i++)
            A[0][ i] = a_big[i];
        for (int i = 0; i < n; i++)
            A[1][ i] = b_big[i];
        for (int i = 0; i < n; i++)
            A[2][ i] = c_big[i];
        // вычисление alpha, beta
        for (int i = 1; i < n; i++)
        {
            alpha[i + 1] = c_big[i] / (b_big[i] - a_big[i] * alpha[i]);
            beta[i + 1] = (a_big[i] * beta[i] - fi[i]) / (b_big[i] - a_big[i] * alpha[i]);
        }
        // вычисление y[i]
        for (int i = n - 1; i > 0; i--)
        {
            y[i] = alpha[i + 1] * y[i + 1] + beta[i + 1];
        }


        y_jacobi = Jacobi(ai, fi, qi);
//        печать метода Якоби
        System.out.println("         y    " + "   y_jacobi   " +"  abs");
        for(int i=0;i<=n;i++){
            System.out.printf("%.1f  %.7f  %.7f  %.7f",i * h,y[i], y_jacobi[i], Math.abs(y[i] - y_jacobi[i]));
            System.out.println();
        }
        System.out.println("Количество итераций: "+" "+k);
        System.out.println();

        y_zeid = Zeid(ai, fi, qi);
        //печать метода Зейделя
        System.out.println("         y    " + "   y_zeid   " +"  abs");
        for(int i=0;i<=n;i++){
            System.out.printf("%.1f  %.7f  %.7f  %.7f",i * h,y[i], y_zeid[i], Math.abs(y[i] - y_zeid[i]));
            System.out.println();
        }
        System.out.println("Количество итераций: "+" "+k);

        y_relax = Relax(ai, fi, qi);
        //печать метода релаксации
        System.out.println("Вычисления для омега={0}"+" "+best_omega);
        System.out.println("         y    " + "   y_relax   " +"  abs");
        for (int i = 0; i <= n; i++) {
            System.out.printf("%.1f  %.7f  %.7f  %.7f", i * h, y[i], y_relax[i], Math.abs(y[i] - y_relax[i]));
            System.out.println();
        }
           System.out.println("Количество итераций: "+" "+k);



    }
    public static double[] Jacobi(double[] a,double[] f,double[] g) {
        double[] y = new double[n + 1];
        double[] y_pre = new double[n + 1];

        double r;
        k = 0;
        for (int i = 0; i < n; i++)
            y[i] = 0;
//        y[i] = Ux(i * h);
//         y[i] = Math.sin(i*h);
        do
        {
            for (int i = 0; i < n; i++)
                y_pre[i] = y[i];
            r = -1;
            for (int i = 1; i < n; i++)
                y[i] = (f[i] + a[i + 1] * y_pre[i + 1] + a[i] * y_pre[i - 1]) / (a[i] + a[i + 1] + g[i] * h * h);
            k++;
            for (int i = 1; i < n; i++)
                if (Math.abs(-a[i + 1] * y[i + 1] + (a[i + 1] + a[i] + g[i] * h * h) * y[i] - a[i] * y[i - 1] - f[i]) > r)
                    r = Math.abs(-a[i + 1] * y[i + 1] + (a[i + 1] + a[i] + g[i] * h * h) * y[i] - a[i] * y[i - 1] - f[i]);
        } while (r > eps);
        return y;
    }

    public static double[] Relax(double[] a, double[] f,double[] g) {
        double[] y = new double[n + 1];
        double[] y_pre = new double[n + 1];
        double r = 0;

        for(double omega = 1.05; omega < 2; omega += 0.05)
        {
            for (int i = 0; i <= n; i++)
                y[i] = 0;
            k = 0;
            do {
                for (int i = 0; i <= n; i++)
                    y_pre[i] = y[i];
                r = -1;
                for (int i = 1; i < n; i++)
                    y[i] = (f[i] + a[i + 1] * y_pre[i + 1] + a[i] * y[i - 1]) / (a[i + 1] + a[i] + g[i] * h * h) * omega + (1 - omega) * y_pre[i];
                for (int i = 1; i < n; i++)
                    if (Math.abs(-a[i + 1] * y[i + 1] + (a[i + 1] + a[i] + g[i] * h * h) * y[i] - a[i] * y[i - 1] - f[i]) > r)
                        r = Math.abs(-a[i + 1] * y[i + 1] + (a[i + 1] + a[i] + g[i] * h * h) * y[i] - a[i] * y[i - 1] - f[i]);
                k++;
            } while (r > eps);

            if(k<best_k) {
                best_omega = omega;
                best_k = k;
            }
            System.out.printf("%.2f  %d",omega,k);
//            System.out.printf("%d  ",k);
            System.out.println();
        }
        for(int i=0;i<=n;i++)
            y[i]=0;
        k=0;
        do {
            for (int i = 0; i <= n; i++)
                y_pre[i] = y[i];
            r = -1;
            for (int i = 1; i < n; i++)
                y[i] = (f[i] + a[i + 1] * y_pre[i + 1] + a[i] * y[i - 1]) / (a[i + 1] + a[i] + g[i] * h * h) * best_omega + (1 - best_omega) * y_pre[i];
            k++;
            for (int i = 1; i < n; i++)
                if (Math.abs(-a[i + 1] * y[i + 1] + (a[i + 1] + a[i] + g[i] * h * h) * y[i] - a[i] * y[i - 1] - f[i]) > r)
                    r = Math.abs(-a[i + 1] * y[i + 1] + (a[i + 1] + a[i] + g[i] * h * h) * y[i] - a[i] * y[i - 1] - f[i]);
        }
        while (r > eps);
        return y;
    }

    public static double[] Zeid(double[] a, double[] f,double[] g) {
        double[] y = new double[n + 1];
        double[] y_pre = new double[n + 1];
        double r = 0;
        double omega=1;
            for (int i = 0; i <= n; i++)
                y[i] = 0;
//              y[i] = Math.sin(i*h);
//              y[i] = Ux(i * h);
            k = 0;
            do {
                for (int i = 0; i <= n; i++)
                    y_pre[i] = y[i];
                r = -1;
                for (int i = 1; i < n; i++)
                    y[i] = (f[i] + a[i + 1] * y_pre[i + 1] + a[i] * y[i - 1]) / (a[i + 1] + a[i] + g[i] * h * h) * omega + (1 - omega) * y_pre[i];
                for (int i = 1; i < n; i++)
                    if (Math.abs(-a[i + 1] * y[i + 1] + (a[i + 1] + a[i] + g[i] * h * h) * y[i] - a[i] * y[i - 1] - f[i]) > r)
                        r = Math.abs(-a[i + 1] * y[i + 1] + (a[i + 1] + a[i] + g[i] * h * h) * y[i] - a[i] * y[i - 1] - f[i]);
                k++;
            } while (r >eps);

        return y;
    }
}

