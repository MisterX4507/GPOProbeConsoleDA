using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading;
using System.Threading.Tasks;
using System.Web.UI.DataVisualization.Charting;
using System.IO;

namespace ПроектГПОКонсоль
{
    delegate double TFunc(double[] x);

    class Program
    {
        public static double[] Aver;
        public static Random rnd11 = new Random();
        static void Main(string[] args)
        {
            int SearchAgents_no, Max_iteration, dim; string Function_name;
            double[] lb; double[] ub; TFunc fobj; double lb1, ub1; double Best_score = 100;
            double[] Best_pos;
            SearchAgents_no = 30; // Number of search agents
            Function_name = "F1"; // Name of the test function that can be from F1 to F13 (Table 1,2,3 in the paper)
            Max_iteration = 500; // Maximum number of iterations
            int Max_progon = 1; Aver = new double[Max_iteration];
            //Load details of the selected benchmark function
            Get_Functions_details(Function_name, out lb1, out ub1, out dim, out fobj);
            for (int k = 0; k < Max_progon; k++)
            {
                DA(SearchAgents_no, Max_iteration, out lb, out ub, dim, fobj, lb1, ub1, out Best_pos, ref Best_score);
                Console.WriteLine("The best optimal value of the objective funciton found by DA is : {0}", Best_score);
            }
            using (StreamWriter sw = new StreamWriter(@"C:\Users\User\Desktop\ll.txt", false))
            {
                for (int k = 0; k < Max_iteration; k++)
                {
                    Aver[k] = Aver[k] / (double)Max_progon;
                    sw.WriteLine(Aver[k]);
                }
            }

            //Console.WriteLine("The best solution obtained by DA is : ");
            //for (int k = 0; k < Best_pos.Length; k++)
            //{
            //    Console.WriteLine(Best_pos[k]);
            //}
            //Console.WriteLine("The best optimal value of the objective funciton found by DA is : {0}", Best_score);
        }
        static void Get_Functions_details(string F, out double lb1, out double ub1, out int dim, out TFunc fobj)
        {
            fobj = F1;
            lb1 = -100;
            ub1 = 100;
            dim = 10;
            switch (F)
            {
                case "F1":
                    break;
                case "F2":
                    fobj = F2;
                    lb1 = -10;
                    ub1 = 10;
                    dim = 10;
                    break;
                case "F3":
                    fobj = F3;
                    lb1 = -100;
                    ub1 = 100;
                    dim = 10;
                    break;
                case "F4":
                    fobj = F4;
                    lb1 = -100;
                    ub1 = 100;
                    dim = 10;
                    break;
                case "F5":
                    fobj = F5;
                    lb1 = -30;
                    ub1 = 30;
                    dim = 10;
                    break;
                case "F6":
                    fobj = F6;
                    lb1 = -100;
                    ub1 = 100;
                    dim = 10;
                    break;
                case "F7":
                    fobj = F7;
                    lb1 = -1.28;
                    ub1 = 1.28;
                    dim = 10;
                    break;
                case "F8":
                    fobj = F8;
                    lb1 = -500;
                    ub1 = 500;
                    dim = 10;
                    break;
                case "F9":
                    fobj = F9;
                    lb1 = -5.12;
                    ub1 = 5.12;
                    dim = 10;
                    break;
                case "F10":
                    fobj = F10;
                    lb1 = -32;
                    ub1 = 32;
                    dim = 10;
                    break;
                case "F11":
                    fobj = F11;
                    lb1 = -600;
                    ub1 = 600;
                    dim = 10;
                    break;
                case "F12":
                    fobj = F12;
                    lb1 = -50;
                    ub1 = 50;
                    dim = 10;
                    break;
                case "F13":
                    fobj = F13;
                    lb1 = -50;
                    ub1 = 50;
                    dim = 10;
                    break;
            }
        }
        static double F1(double[] x)
        {
            double sum = 0; foreach (double elem in x)
            {
                sum = sum + Math.Pow(elem, 2);
            }
            return sum;
        }
        static double F2(double[] x)
        {
            double sum = 0, prod = 1; foreach (double elem in x)
            {
                sum = sum + Math.Abs(elem);
                prod = prod * Math.Abs(elem);
            }
            return (sum + prod);
        }
        static double F3(double[] x)
        {
            int i, j, dim = x.Length; double y = 0;
            for (i = 0; i < dim; i++)
            {
                double sum = 0;
                for (j = 0; j <= i; j++)
                {
                    sum = sum + x[j];
                }
                y = y + Math.Pow(sum, 2);
            }
            return y;
        }
        static double F4(double[] x)
        {
            double max = 0;
            foreach (double elem in x)
            {
                if (max < Math.Abs(elem)) max = Math.Abs(elem);
            }
            return max;
        }
        static double F5(double[] x)
        {
            int dim = x.Length, i; double sum = 0;
            for (i = 0; i <= (dim - 2); i++)
            {
                sum = sum + (100 * Math.Pow((x[i + 1] - Math.Pow(x[i], 2)), 2) + Math.Pow((x[i] - 1), 2));
            }
            return sum;
        }
        static double F6(double[] x)
        {
            double sum = 0; foreach (double elem in x)
            {
                sum = sum + Math.Pow((elem + 0.5), 2);
            }
            return sum;
        }
        static double F7(double[] x)
        {
            double sum = rnd11.NextDouble(); int i, dim = x.Length;
            for (i = 0; i < dim; i++)
            {
                sum = sum + (i + 1) * Math.Pow(x[i], 4);
            }
            return sum;
        }
        static double F8(double[] x)
        {
            double sum = 0; foreach (double elem in x)
            {
                sum = sum + (-1 * elem * Math.Sin(Math.Sqrt(Math.Abs(elem))));
            }
            return sum;
        }
        static double F9(double[] x)
        {
            double sum = 0; foreach (double elem in x)
            {
                sum = sum + (Math.Pow(elem, 2) - 10 * Math.Cos(2 * Math.PI * elem) + 10);
            }
            return sum;
        }
        static double F10(double[] x)
        {
            double sum1 = 0, sum2 = 0, sum;
            foreach (double elem in x)
            {
                sum1 = sum1 + Math.Pow(elem, 2);
                sum2 = sum2 + Math.Cos(2 * Math.PI * elem);
            }
            sum1 = sum1 / x.Length; sum2 = sum2 / x.Length;
            sum = -20 * Math.Exp(-0.2 * Math.Sqrt(sum1)) - Math.Exp(sum2) + 20 + Math.Exp(1);
            return sum;
        }
        static double F11(double[] x)
        {
            double sum = 0, mply = 1; int i, dim = x.Length;
            for (i = 0; i < dim; i++)
            {
                sum = sum + Math.Pow(x[i], 2);
                mply = mply * Math.Cos(x[i] / Math.Sqrt(i + 1));
            }
            double y = sum / 4000 - mply + 1;
            return y;
        }
        static double F12(double[] x)
        {
            int dim = x.Length;
            double[] z = new double[dim];
            for (int i = 0; i < dim; i++)
            {
                z[i] = 1 + (x[i] + 1) / 4.0;
            }
            double sum1 = 0, sum2 = 0;
            for (int i = 0; i < (dim - 1); i++)
            {
                sum1 = sum1 + Math.Pow(z[i] - 1, 2) * (1 + 10 * Math.Pow(Math.Sin(Math.PI * z[i + 1]), 2));
            }
            foreach (double elem in x)
            {
                sum2 = sum2 + uFunc(elem, 10, 100, 4);
            }
            double y = (Math.PI / dim) * (10 * Math.Sin(Math.PI * z[0]) + sum1 + Math.Pow(z[dim - 1] - 1, 2)) + sum2;
            return y;
        }
        static double uFunc(double x, int a, int k, int m)
        {
            if (x > a) return (k * Math.Pow(x - a, m));
            if (x < -a) return (k * Math.Pow(-x - a, m));
            return 0;
        }
        static double F13(double[] x)
        {
            int dim = x.Length; double sum1 = 0, sum2 = 0;
            foreach (double elem in x)
            {
                sum1 = sum1 + Math.Pow(elem - 1, 2) * (1 + Math.Pow(Math.Sin(3 * Math.PI * elem + 1), 2));
                sum2 = sum2 + uFunc(elem, 5, 100, 4);
            }
            double y = 0.1 * (Math.Pow(Math.Sin(3 * Math.PI * x[0]), 2) + sum1 + Math.Pow(x[dim - 1] - 1, 2) * (1 + Math.Pow(Math.Sin(2 * Math.PI * x[dim - 1]), 2))) + sum2;
            return y;
        }
        static void initialization(int SearchAgents_no, int dim, double[] lb, double[] ub, ref double[][] dX)
        {
            int i, j;
            for (i = 0; i < SearchAgents_no; i++)
            {
                for (j = 0; j < dim; j++)
                {
                    dX[i][j] = rnd11.NextDouble() * (ub[j] - lb[j]) + lb[j];
                }
            }
        }
        static void DA(int SearchAgents_no, int Max_iteration, out double[] lb, out double[] ub, int dim, TFunc fobj, double lb1, double ub1, out double[] Best_p, ref double Best_s)
        {
            int i, j; lb = new double[dim]; ub = new double[dim]; Best_p = new double[dim];
            for (i = 0; i < dim; i++)
            {
                lb[i] = lb1; ub[i] = ub1;
            }
            double[] r = new double[dim], Delta_max = new double[dim];
            for (i = 0; i < dim; i++)
            {
                r[i] = (ub[i] - lb[i]) / 10.0;
                Delta_max[i] = (ub[i] - lb[i]) / 10.0; // The initial radius of gragonflies' neighbourhoods
            }
            double Food_fitness, Enemy_fitness;
            double[] Food_pos = new double[dim], Enemy_pos = new double[dim];
            Food_fitness = Double.PositiveInfinity; Enemy_fitness = Double.NegativeInfinity;
            double[][] X = new double[SearchAgents_no][]; double[][] DeltaX = new double[SearchAgents_no][];
            for (i = 0; i < SearchAgents_no; i++)
            {
                X[i] = new double[dim];
                DeltaX[i] = new double[dim];
            }
            initialization(SearchAgents_no, dim, lb, ub, ref X);
            double[] Fitness = new double[SearchAgents_no];
            initialization(SearchAgents_no, dim, lb, ub, ref DeltaX);
            int iter; double my_c, w;
            for (iter = 1; iter <= Max_iteration; iter++)
            {
                for (i = 0; i < dim; i++)
                {
                    r[i] = (ub[i] - lb[i]) / 4.0 + ((ub[i] - lb[i]) * ((double)iter / (double)Max_iteration) * 2);
                }
                w = 0.9 - iter * (0.5 / (double)Max_iteration);
                my_c = 0.1 - iter * (0.1 / ((double)Max_iteration / 2.0));
                //Console.WriteLine("w={0}   my_c={1}", w, my_c); Console.Read();
                if (my_c < 0) my_c = 0;
                double s, a, c, f, e;
                s = 2 * rnd11.NextDouble() * my_c; // Seperation weight                
                a = 2 * rnd11.NextDouble() * my_c; // Alignment weight                
                c = 2 * rnd11.NextDouble() * my_c; // Cohesion weight         
                f = 2 * rnd11.NextDouble();        // Food attraction weight
                e = my_c;                         // Enemy distraction weight
                for (i = 0; i < SearchAgents_no; i++) //Calculate all the objective values
                {
                    Fitness[i] = fobj(X[i]);
                    //Console.Write("Fitness[{0}][{1}]={2}", iter, i, Fitness[i]); Console.ReadLine();
                    if (Fitness[i] < Food_fitness)
                    {
                        Food_fitness = Fitness[i];
                        //Console.WriteLine("{0}: {1}", iter, Food_fitness); Console.Read();
                        for (j = 0; j < dim; j++)
                        {
                            Food_pos[j] = X[i][j];
                            //Console.WriteLine("Food_pos[{0}]={1}", j, Food_pos[j]); 
                        }
                    }
                    bool q = true;
                    if (Fitness[i] > Enemy_fitness)
                    {
                        for (j = 0; j < dim; j++)
                        {
                            if ((X[i][j] > ub[j]) || (X[i][j] < lb[j])) { q = false; break; }
                        }
                        if (q == true)
                        {
                            Enemy_fitness = Fitness[i];
                            for (j = 0; j < dim; j++)
                            {
                                Enemy_pos[j] = X[i][j];
                            }
                        }
                    }
                }
                for (i = 0; i < SearchAgents_no; i++)
                {
                    double Dist2Enemy = 0;
                    double Dist2Food = 0;
                    int index = 0; int neighbours_no = 0;
                    double[][] Neighbours_DeltaX = new double[SearchAgents_no][]; //clear?
                    double[][] Neighbours_X = new double[SearchAgents_no][];
                    //find the neighbouring solutions
                    for (j = 0; j < SearchAgents_no; j++)
                    {
                        distance(X[i], X[j], ref Dist2Enemy);
                        //Console.WriteLine("{0} и {1}: {2}", i, j, Dist2Enemy);
                        bool p = true; for (int k = 0; k < dim; k++)
                        {
                            if ((Dist2Enemy > r[k]) || (Dist2Enemy == 0))
                            {
                                p = false; break;
                            }
                        }
                        if (p == true)
                        {
                            //if ((i == 29) && (j == 27))
                            //{
                            //    for (int k = 0; k < dim; k++)
                            //    {
                            //        Console.Write("{0} ", X[i][k]);
                            //    }
                            //    Console.WriteLine();
                            //    for (int k = 0; k < dim; k++)
                            //    {
                            //       Console.Write("{0} ", X[j][k]);
                            //   }
                            //    Console.WriteLine();
                            //}
                            //Console.WriteLine("Растояние между {0} и {1} = {2}, радиус = {3}", i, j, Dist2Enemy, r[0]);
                            index = index + 1; neighbours_no = neighbours_no + 1;
                            Neighbours_DeltaX[index - 1] = new double[dim];
                            Neighbours_X[index - 1] = new double[dim];
                            for (int k = 0; k < dim; k++)
                            {
                                //Console.WriteLine("{4} стрекоза: DeltaX[{0}][{1}]={2} и X[{0}][{1}]={3}", j, k, DeltaX[j][k], X[j][k], i);
                                Neighbours_DeltaX[index - 1][k] = DeltaX[j][k];  //??????????????????????????
                                Neighbours_X[index - 1][k] = X[j][k];
                                //Console.WriteLine("Neighbours_DeltaX[{0}][{1}]={2}; Neighbours_X[{0}][{1}]={3}", index - 1, k, Neighbours_DeltaX[index - 1][k], Neighbours_X[index - 1][k]);
                            }
                        }
                    }
                    //Console.WriteLine("{0}: Число соседей = {1}", i, neighbours_no);
                    //for (int yyy = 0; yyy < neighbours_no; yyy++)
                    //{
                    //    for (int k = 0; k < dim; k++)
                    //    {
                    //        Console.WriteLine("Neighbours_DeltaX[{0}][{1}]={2}; Neighbours_X[{0}][{1}]={3}", yyy, k, Neighbours_DeltaX[yyy][k], Neighbours_X[yyy][k]);
                    //    }
                    //}                        
                    // Seperation Eq. (3.1)
                    double[] S = new double[dim];
                    if (neighbours_no >= 1)
                    {
                        for (int k = 0; k < neighbours_no; k++)
                        {
                            for (int ka = 0; ka < dim; ka++)
                            {
                                S[ka] = S[ka] + (Neighbours_X[k][ka] - X[i][ka]);  // mistakable?
                            }
                        }
                        for (int k = 0; k < dim; k++)
                        {
                            S[k] = -1 * S[k];
                        }
                    }
                    else
                    {
                        Array.Clear(S, 0, dim);
                    }
                    //for (int k = 0; k < dim; k++)
                    //{
                    //    Console.Write("S[{0}][{1}]={2} ", i, k, S[k]);
                    //}
                    //Console.ReadLine();
                    // Alignment Eq. (3.2)
                    double[] A = new double[dim];
                    if (neighbours_no > 1)
                    {
                        for (int k1 = 0; k1 < dim; k1++)
                        {
                            for (int k2 = 0; k2 < neighbours_no; k2++)
                            {
                                A[k1] = A[k1] + Neighbours_DeltaX[k2][k1];
                            }
                            A[k1] = A[k1] / (double)neighbours_no;
                        }
                    }
                    else
                    {
                        for (int k1 = 0; k1 < dim; k1++)
                        {
                            A[k1] = DeltaX[i][k1];
                        }
                    }
                    // Cohesion Eq. (3.3)
                    double[] C_temp = new double[dim];
                    if (neighbours_no >= 1)
                    {
                        for (int k1 = 0; k1 < dim; k1++)
                        {
                            for (int k2 = 0; k2 < neighbours_no; k2++)
                            {
                                C_temp[k1] = C_temp[k1] + Neighbours_X[k2][k1];
                            }
                            C_temp[k1] = C_temp[k1] / (double)neighbours_no;
                        }
                    }
                    else
                    {
                        for (int k1 = 0; k1 < dim; k1++)
                        {
                            C_temp[k1] = X[i][k1];
                        }
                    }
                    double[] C = new double[dim];
                    for (int k = 0; k < dim; k++)
                    {
                        C[k] = C_temp[k] - X[i][k];
                    }
                    for (int k = 0; k < dim; k++)
                    {
                        Console.Write("C[{0}][{1}]={2} ", i, k, C[k]);
                    }
                    Console.WriteLine();
                    // Attraction to food Eq. (3.4)
                    distance(X[i], Food_pos, ref Dist2Food);
                    //Console.WriteLine("Еда {0}: {1}", i, Dist2Food);
                    double[] F = new double[dim];
                    bool q1 = true; for (int k = 0; k < dim; k++)
                    {
                        if (Dist2Food > r[k]) { q1 = false; break; }
                    }
                    if (q1 == true)
                    {
                        for (int k = 0; k < dim; k++)
                        {
                            F[k] = Food_pos[k] - X[i][k];
                        }
                    }
                    else
                    {
                        Array.Clear(F, 0, dim);
                    }
                    // Distraction from enemy Eq. (3.5)
                    distance(X[i], Enemy_pos, ref Dist2Enemy);
                    //Console.WriteLine("Враг {0}: {1}", i, Dist2Enemy);
                    double[] Enemy = new double[dim];
                    bool q2 = true; for (int k = 0; k < dim; k++)
                    {
                        if (Dist2Enemy > r[k]) { q2 = false; break; }
                    }
                    if (q2 == true)
                    {
                        for (int k = 0; k < dim; k++)
                        {
                            Enemy[k] = Enemy_pos[k] + X[i][k];
                        }
                    }
                    else
                    {
                        Array.Clear(Enemy, 0, dim);
                    }
                    //if (i == 2)
                    //{
                    //Console.WriteLine("{0}, {1}:", iter, i);
                    //for (int tt = 0; tt<dim; tt++)
                    //{
                    //Console.Write("S[{0}]={1} ", tt, S[tt]);
                    //}
                    //Console.WriteLine();
                    //for (int tt = 0; tt < dim; tt++)
                    //{
                    //    Console.Write("A[{0}]={1} ", tt, A[tt]);
                    //}
                    //Console.WriteLine();
                    //for (int tt = 0; tt < dim; tt++)
                    //{
                    //    Console.Write("C[{0}]={1} ", tt, C[tt]);
                    //}
                    //Console.WriteLine();
                    //for (int tt = 0; tt < dim; tt++)
                    //{
                    //    Console.Write("F[{0}]={1} ", tt, F[tt]);
                    //}
                    //Console.WriteLine();
                    //for (int tt = 0; tt < dim; tt++)
                    //{
                    //    Console.Write("Enemy[{0}]={1} ", tt, Enemy[tt]);
                    //}
                    //Console.WriteLine();
                    //}
                    for (int tt = 0; tt < dim; ++tt)
                    {
                        if (X[i][tt] > ub[tt])
                        {
                            X[i][tt] = lb[tt];
                            DeltaX[i][tt] = rnd11.NextDouble();
                        }
                        if (X[i][tt] < lb[tt])
                        {
                            X[i][tt] = ub[tt];
                            DeltaX[i][tt] = rnd11.NextDouble();
                        }
                    }
                    q1 = true; for (int k = 0; k < dim; k++)
                    {
                        if (Dist2Food > r[k]) { q1 = false; break; }
                    }
                    if (q1 == false)
                    {
                        if (neighbours_no >= 1)
                        {
                            double r1, r2, r3;
                            r1 = rnd11.NextDouble();
                            r2 = rnd11.NextDouble();
                            r3 = rnd11.NextDouble();
                            for (j = 0; j < dim; j++)
                            {
                                DeltaX[i][j] = w * DeltaX[i][j] + r1 * S[j] + r2 * A[j] + r3 * C[j];
                                if (DeltaX[i][j] > Delta_max[j])
                                {
                                    DeltaX[i][j] = Delta_max[j];
                                }
                                if (DeltaX[i][j] < (-1 * Delta_max[j]))
                                {
                                    DeltaX[i][j] = (-1 * Delta_max[j]);
                                }
                                X[i][j] = X[i][j] + DeltaX[i][j];
                            }
                        }
                        else
                        {
                            // Eq. (3.8)
                            double[] Ly = new double[dim]; Ly = Levy(dim);
                            for (j = 0; j < dim; j++)
                            {
                                X[i][j] = X[i][j] + Ly[j] * X[i][j];
                                DeltaX[i][j] = 0;
                            }
                        }
                    }
                    else
                    {
                        for (j = 0; j < dim; j++)
                        {
                            // Eq. (3.6)
                            DeltaX[i][j] = (s * S[j] + a * A[j] + c * C[j] + f * F[j] + e * Enemy[j]) + w * DeltaX[i][j];
                            if (DeltaX[i][j] > Delta_max[j])
                            {
                                DeltaX[i][j] = Delta_max[j];
                            }
                            if (DeltaX[i][j] < (-1 * Delta_max[j]))
                            {
                                DeltaX[i][j] = -1 * Delta_max[j];
                            }
                            X[i][j] = X[i][j] + DeltaX[i][j];
                        }
                    }
                    //if (i == 1) { Console.WriteLine(iter); }
                    for (j = 0; j < dim; j++)
                    {
                        if (X[i][j] > ub[j])
                        {
                            X[i][j] = ub[j];
                        }
                        if (X[i][j] < lb[j])
                        {
                            X[i][j] = lb[j];
                        }
                        //if (i == 1)
                        //{ Console.Write("{0} ", DeltaX[i][j]); }
                    }
                    //if (i == 1)
                    // Console.ReadLine();
                }
                Best_s = Food_fitness;
                Aver[iter - 1] = Aver[iter - 1] + Food_fitness;
                Best_p = Food_pos;
            }
        }
        static void distance(double[] A, double[] B, ref double Dist)
        {
            Dist = 0;
            for (int i = 0; i < A.Length; i++)
            {
                Dist = Dist + Math.Pow(A[i] - B[i], 2);
            }
            Dist = Math.Sqrt(Dist);
        }
        static double[] Levy(int d)
        {
            double beta = 1.5; double sigma, r1, r2, r3, r4;
            Chart Chart1 = new Chart();
            // Eq. (3.10)
            sigma = Math.Pow((Chart1.DataManipulator.Statistics.GammaFunction(1 + beta) * Math.Sin(Math.PI * beta / 2.0) / (Chart1.DataManipulator.Statistics.GammaFunction((1 + beta) / 2.0) * beta * Math.Pow(2, ((beta - 1) / 2.0)))), (1.0 / beta));
            double[] step = new double[d];
            for (int i = 0; i < d; i++)
            {
                r1 = rnd11.NextDouble();
                r2 = rnd11.NextDouble();
                r3 = rnd11.NextDouble();
                r4 = rnd11.NextDouble();
                step[i] = ((randn(r1, r3) * sigma) / Math.Pow(Math.Abs(randn(r2, r4)), (1.0 / beta))) * 0.01;
            }
            return step;
        }
        static double randn(double r, double p)
        {
            double n;
            n = Math.Sqrt(-1 * 2 * Math.Log(r)) * Math.Cos(2 * Math.PI * p);
            return n;
        }
    }
}