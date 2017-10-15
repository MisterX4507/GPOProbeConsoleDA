﻿using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace ПроектГПОКонсоль
{
    delegate double TFunc(double[] x);

    class Program
    {
        static void Main(string[] args)
        {
            int SearchAgents_no, Max_iteration, dim; string Function_name;
            double[] lb; double[] ub; TFunc fobj; double lb1, ub1;
            SearchAgents_no = 40; // Number of search agents
            Function_name = "F1"; // Name of the test function that can be from F1 to F23 (Table 1,2,3 in the paper)
            Max_iteration = 500; // Maximum numbef of iterations
            //Load details of the selected benchmark function
            Get_Functions_details(Function_name, out lb1, out ub1, out dim, out fobj);
            DA(SearchAgents_no, Max_iteration, out lb, out ub, dim, fobj, lb1, ub1);
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
            for (i = 0; i < dim; ++i)
            {
                double sum = 0;
                for (j = 0; j <= i; ++j)
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
            for (i = 0; i <= (dim - 2); ++i)
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
            Random rnd = new Random();
            double sum = rnd.NextDouble(); int i, dim = x.Length;
            for (i = 0; i < dim; ++i)
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
            for (i = 0; i < dim; ++i)
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
            for (int i = 0; i < dim; ++i)
            {
                z[i] = 1 + (x[i] + 1) / 4.0;
            }
            double sum1 = 0, sum2 = 0;
            for (int i = 0; i < (dim - 1); ++i)
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
        static void initialization(int SearchAgents_no, int dim, double[] lb, double[] ub, ref double[][] X)
        {
            int i, j; Random rnd = new Random();
            for (i = 0; i < SearchAgents_no; ++i)
            {
                for (j = 0; j < dim; ++j)
                {
                    X[i][j] = rnd.NextDouble() * (ub[i] - lb[i]) + lb[i];
                }
            }
        }
        static void DA(int SearchAgents_no, int Max_iteration, out double[] lb, out double[] ub, int dim, TFunc fobj, double lb1, double ub1)
        {
            int i, j; lb = new double [dim]; ub = new double [dim];
            for (i = 0; i < dim; ++i)
            {
                lb[i] = lb1; ub[i] = ub1;
            }
            double[] r=new double [dim], Delta_max=new double [dim];
            for (i = 0; i < dim; ++i)
            {
                Delta_max[i] = (ub[i] - lb[i]) / 10.0; // The initial radius of gragonflies' neighbourhoods
            }
            double Food_fitness, Enemy_fitness;
            double[] Food_pos = new double[dim], Enemy_pos = new double[dim];
            Food_fitness = Double.PositiveInfinity; Enemy_fitness = Double.NegativeInfinity;
            double[][] X = new double[SearchAgents_no][]; double[][] DeltaX = new double[SearchAgents_no][];
            for (i = 0; i < SearchAgents_no; ++i)
            {
                X[i] = new double[dim];
                DeltaX[i] = new double[dim];
            }
                initialization(SearchAgents_no, dim, lb, ub, ref X);
            double[] Fitness = new double[SearchAgents_no];
            initialization(SearchAgents_no, dim, lb, ub, ref DeltaX);
            int iter; double my_c, w;
            for (iter = 1; iter <= Max_iteration; ++i)
            {
                for (i = 0; i < dim; ++i)
                {
                    r[i] = (ub[i] - lb[i]) / 4.0 + ((ub[i]-lb[i])*((double)iter/(double)Max_iteration)*2);
                }
                w = 0.9 - iter * (0.5 / (double)Max_iteration);
                my_c = 0.1 - iter * (0.1 / ((double)Max_iteration / 2.0));
                if (my_c < 0) my_c = 0;
                double s, a, c, f, e; Random rnd = new Random();
                s = 2 * rnd.NextDouble() * my_c; // Seperation weight
                a = 2 * rnd.NextDouble() * my_c; // Alignment weight
                c = 2 * rnd.NextDouble() * my_c; // Cohesion weight
                f = 2 * rnd.NextDouble();        // Food attraction weight
                e = my_c;                        // Enemy distraction weight
                for (i = 0; i < SearchAgents_no; ++i) //Calculate all the objective values
                {                    
                    Fitness[i] = fobj(X[i]);
                    if (Fitness[i] < Food_fitness)
                    {
                        Food_fitness = Fitness[i];
                        Food_pos = X[i];
                    }
                    bool q=true;
                    if (Fitness[i] > Enemy_fitness)
                    {
                        for (j = 0; j < dim; ++j)
                        {
                            if ((X[i][j] > ub[j]) || (X[i][j] < lb[j])) { q = false; break; }
                        }
                        if (q)
                        {
                            Enemy_fitness = Fitness[i];
                            Enemy_pos = X[i];
                        }
                    }
                }
                for (i = 0; i < SearchAgents_no; ++i)
                {
                    double[] Dist2Enemy=new double [dim];
                    double[] Dist2Food = new double [dim];
                    int index = 0; int neighbours_no = 0;
                    double[][] Neighbours_DeltaX = new double[SearchAgents_no][]; //clear?
                    double[][] Neighbours_X = new double[SearchAgents_no][];
                    //find the neighbouring solutions
                    for (j = 0; j < SearchAgents_no; ++j)
                    {
                        distance(X[i], X[j], ref Dist2Enemy);
                        bool q = true; for (int k = 0; k < dim; ++k)
                        {
                            if ((Dist2Enemy[k] > r[k]) || (Dist2Enemy[k] == 0)) { q = false; break; }
                        }
                        if (q)
                        {
                            index = index + 1; neighbours_no = neighbours_no + 1;
                            Neighbours_DeltaX[index - 1] = new double[dim];
                            Neighbours_X[index - 1] = new double[dim];
                            Neighbours_DeltaX[index-1] = DeltaX[j];
                            Neighbours_X[index-1] = X[j];
                        }
                    }
                    // Seperation Eq. (3.1)
                    double[] S = new double[dim];
                    if (neighbours_no > 1)
                    {
                        for (int k = 0; k < neighbours_no; ++k)
                        {
                            for (int ka = 0; ka < dim; ++ka)
                            {
                                S[ka] = S[ka] + (Neighbours_X[k][ka] - X[i][ka]);  // mistakable?
                            }
                        }
                        for (int k = 0; k < S.Length; ++k)
                        {
                            S[k] = -1 * S[k];
                        }
                    }
                    else
                    {
                        Array.Clear(S, 0, dim);
                    }
                    // Alignment Eq. (3.2)
                    double[] A = new double[dim];
                    if (neighbours_no > 1)
                    {
                        for (int k1 = 0; k1 < dim; ++k1)
                        {
                            for (int k2 = 0; k2 < index; ++k2)
                            {
                                A[k1] = A[k1] + Neighbours_DeltaX[k2][k1];
                            }
                            A[k1] = A[k1] / (double)neighbours_no;
                        }
                    }
                    else
                    {
                        A = DeltaX[i];
                    }
                    // Cohesion Eq. (3.3)
                    double[] C_temp = new double[dim];
                    if (neighbours_no > 1)
                    {
                        for (int k1 = 0; k1 < dim; ++k1)
                        {
                            for (int k2 = 0; k2 < index; ++k2)
                            {
                                C_temp[k1] = C_temp[k1] + Neighbours_X[k2][k1];
                            }
                            C_temp[k1] = C_temp[k1] / (double)neighbours_no;
                        }
                    }
                    else
                    {
                        C_temp = X[i];
                    }
                    double[] C = new double[dim];
                    for (int k = 0; k < dim; ++k)
                    {
                        C[k] = C_temp[k] - X[i][k];
                    }
                    // Attraction to food Eq. (3.4)
                    distance(X[i], Food_pos, ref Dist2Food);
                    double[] F = new double[dim];
                    bool q1 = true; for (int k = 0; k < dim; ++k)
                    {
                        if (Dist2Food[k] > r[k]) { q1 = false; break; }
                    }
                    if (q1)
                    {
                        for (int k = 0; k < dim; ++k)
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
                    double[] Enemy = new double[dim];
                    q1 = true; for (int k = 0; k < dim; ++k)
                    {
                        if (Dist2Enemy[k] > r[k]) { q1 = false; break; }
                    }
                    if (q1)
                    {
                        for (int k = 0; k < dim; ++k)
                        {
                            Enemy[k] = Enemy_pos[k] + X[i][k];
                        }
                    }
                    else
                    {
                        Array.Clear(Enemy, 0, dim);
                    }
                    Random rand = new Random();
                    for (int tt = 0; tt < dim; ++tt)
                    {
                        if (X[i][tt] > ub[tt])
                        {
                            X[i][tt] = lb[tt];
                            DeltaX[i][tt] = rand.NextDouble();
                        }
                        if (X[i][tt] < lb[tt])
                        {
                            X[i][tt] = ub[tt];
                            DeltaX[i][tt] = rand.NextDouble();
                        }
                    }

                }
            }
        }
        static void distance(double[] A, double[] B, ref double[] Dist)
        {
            for (int i = 0; i < Dist.Length; ++i)
            {
                Dist[i] = Math.Sqrt(Math.Pow(A[i]-B[i],2));
            }
        }
    }
}