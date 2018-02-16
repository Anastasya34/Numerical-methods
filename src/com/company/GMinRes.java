package com.company;

import Jama.Matrix;
import Jama.LUDecomposition;
import Jama.Matrix;
import Jama.QRDecomposition;
import javax.swing.plaf.multi.MultiViewportUI;
import java.util.*;


public class GMinRes {

    /* Дополнительные методы для Jama.Matrix */
    public static void printMatrix(Matrix matrix){
        for (int i = 0; i < matrix.getRowDimension(); i++){
            for (int j = 0; j < matrix.getColumnDimension(); j++){
                System.out.print(matrix.get(i,j)+" ");
            }
            System.out.println();
        }
    }
    public static mVector mul_Matr_on_Vect(Matrix matrix, mVector vector) throws Exception {
        if (matrix.getColumnDimension() != vector.getLength()) {
            throw new Exception("Размерности матрицы и ветора не совпадают.");
        }
        mVector result = new mVector(vector.getLength());
        for (int i = 0; i < matrix.getRowDimension(); i++) {
            result.vector[i] = 0;
            for (int j = 0; j < matrix.getColumnDimension(); j++) {
                result.vector[i] += matrix.get(i,j) * vector.vector[j];
            }
        }
        return result;
    }
    /* Метод получения псевдослучайного целого числа от min до max (включая max)*/
    public static int rnd(int min, int max){
        max -= min;
        return (int) (Math.random() * ++max) + min;
    }
    /* Генерация данных */
    public static Matrix generateMatrix(int n){
        double [][] matrix = new double[n][n];
        for (int i=0; i < n; i++){
            for (int j=0; j < n; j++){
                if (i == j){
                    matrix[i][i] = rnd(900,5000);
                }
                else{
                    matrix[i][j] = rnd(-5,0);
                }
            }

        }
        return new Matrix(matrix);
    }
    public static Matrix generateMatrixBadCond(int n){
        double [][] matrix = new double[n][n];
        for (int i = 0; i < n; i++){
            int k = rnd(1,2);
            for (int j = 0; j < n; j++){
                if (j % k == 0) {
                    if (i == j){
                        matrix[i][i] = rnd(5,50);
                    }
                    else {
                        matrix[i][j] = 0.0;
                    }
                }
                else {
                    if (i == j){
                        matrix[i][i] = rnd(2000,50000);
                    }
                    else {
                        matrix[i][j] = rnd(-5, 0);
                    }

                }
            }

        }
        return new Matrix(matrix);
    }

    public static mVector generateVector(int n) throws Exception {
        mVector vector = new mVector(n);
        for (int i=0; i < n; i++){
            vector.vector[i] = rnd(0,5000);
        }
        return vector;
    }
    /*Реализация*/
    public static ILUResult ILU(Matrix A, int P){
        Matrix L = new Matrix(A.getRowDimension(), A.getRowDimension(),0.0);
        Matrix U = new Matrix(A.getRowDimension(), A.getRowDimension(),0.0);
        int n = A.getColumnDimension();
        //создание портрета матрицы:
        IndexSet portraitMatrixA = new IndexSet(A);
        for (int p = 0; p <= P; p++) {
            for (int i = 0; i < n; i++) {
                //заполнение первой строки U и первого столбца L
                if (portraitMatrixA.containsPair(0, i)) {
                    //U[0][i] = A[0][i];
                    U.set(0, i,
                            A.get(0, i));
                }
                if (portraitMatrixA.containsPair(i, 0)) {
                    //L[i][0] = A[0][i]/U[0][0];
                    L.set(i, 0,
                            A.get(i, 0) / U.get(0, 0));
                }
                for (int j = i; 0 < j && j < n; j++) {
                    double sum = 0;
                    //заполнение U
                    if (portraitMatrixA.containsPair(i, j)) {
                        for (int s = 0; s < i; s++) {
                            sum += L.get(i, s) * U.get(s, j);
                        }
                        //U[i][j] = A[i][j] - ∑L[i][s]*U[s][j] s = 0..i-1
                        U.set(i, j, A.get(i, j) - sum);
                    }
                    //заполнение L
                    if (portraitMatrixA.containsPair(j, i) && j < n) {
                        sum = 0;
                        for (int t = 0; t < i; t++) {
                            sum += L.get(j, t) * U.get(t, i);
                        }
                        //L[j][i] = (A[j][i] - ∑L[j][t]*U[t][i])/ U[i][i] t = 0..i-1
                        L.set(j, i, (A.get(j, i) - sum) / U.get(i, i));
                    }
                }
            }
            Matrix R = A.minus(L.times(U));
            portraitMatrixA.addPortraitMatrix(R);
        }
        return new ILUResult(L,U);
    }

    private static ArnoldiResult arnoldiIteration(Matrix A, mVector v,Matrix M_inv) throws Exception {
        int n = A.getRowDimension(); //Shape of the input matrix
        Matrix  H = new Matrix(n+1,n,0.0);//Creats a zero matrix of shape (n + 1) x n
        Matrix  V = new Matrix(n,n+1,0.0);//Creats a zero matrix of shape m x n
        mVector w ,z;
        //Adds v to the first column of V
        for (int i = 0; i < n; i++) {
            V.set(i,0,
                    v.get(i));
        }
        for (int i = 0; i < n; i++) {
            //z = A * v_i
            z = mul_Matr_on_Vect(A, v);
            //w = M^(-1)z с учетом предобуславливания
            w = mul_Matr_on_Vect(M_inv,z);
            for (int j = 0; j < i + 1; j++) {
                mVector v_j = new mVector(V, j);
                //h[j][i] = (w,v_j)
                double h_ji = w.mul(v_j);
                H.set(j, i,
                        h_ji);
                //w = w - h[j][i]*v_j
                w = w.subtract(v_j.mul(h_ji));
            }
            // H[i+1][i] = ||w||
            H.set(i + 1, i, w.norm());
            //v_j+1 = w/H[i+1][i]
            v = w.div(H.get(i + 1, i));
            //V[:,j + 1] = v_j+1;
            for (int k = 0; k < n ; k++) {
                V.set(k, i+1,
                        v.get(k));
            }
        }
        return new ArnoldiResult(H,V);
    }


    private static mVector gMRESwithILUcond(Matrix A, mVector b, int p) throws Exception {
        int n = A.getColumnDimension();
        double eps = 0.000000001;
        mVector Xk = new mVector(n, 0.0);
        mVector w = new mVector(n, 0.0);
        Matrix H ;// matrix of shape (n + 1) x n
        Matrix V ;// matrix of shape (n + 1) x n
        mVector e = new mVector(n + 1, 0.0);
        e.set(0, 1); //e = [1,0,...,0];
        mVector r,z;
        int counter = 0;
        ILUResult ILU_decomp = ILU(A, p);
        Matrix M = ILU_decomp.getLU();// M = LU - матрица предобуславливателя
        Matrix M_inv = M.inverse();
        //z = b_o - A X_o
        z = b.subtract(mul_Matr_on_Vect(A,Xk));
        //r_o = M^(-1) z - c учетом предобуславлиявания
        r = mul_Matr_on_Vect(M_inv,z);//невязка начального приближения x0.
        while (r.norm() > eps){
            //w = r/||r|| - начальное значение вектора для построения ортогонального базиса , ||w||=1
            for (int i = 0; i < n; i++) {
                double value = r.get(i) / r.norm();
                w.set(i, value);
            }
            ArnoldiResult arnoldiResult = arnoldiIteration(A, w, M_inv);// найти ортогональный базис Km
            H = arnoldiResult.H;// матрица коэффициентов ортогонализации, расширенная (n+1) x n
            V = arnoldiResult.V;//ортонормированный базис подпространства Km + последний вычисленный вектор v_n. n x (n+1)
            //H = QR
            QRDecomposition QR = new QRDecomposition(H);
            Matrix Q = QR.getQ();
            Matrix R = QR.getR();
            //т.к матрица ортогонвльная Q* = Q^(-1)
            Matrix Q_conj = Q.inverse();
            // Q_conj H = R
            //||Hy - βe || = || Q_conj *(Hy - βe) || = || Ry - β Q_conj e|| = ||Ry - g || = 0
            //β = ||r||
            double beta = r.norm();
            //g = β Q_conj e;
            mVector g = mul_Matr_on_Vect(Q_conj.times(beta), e);
            // y = R^(-1)*g
            mVector y = mul_Matr_on_Vect(R.inverse(),g.getVector(0,n-1));
            //Xk+1 = Xk + V*y
            Xk = Xk.add(mul_Matr_on_Vect(V.getMatrix(0, n - 1, 0, n - 1), y));
            //r = b - A*X_k
            r = b.subtract(mul_Matr_on_Vect(A, Xk));
            counter++;
        }
        System.out.println("Количество итераций метода минимальных невязок: " + counter);
        return Xk;
    }

    public static double cond(Matrix A){
        double normA = A.normInf();
        Matrix Ainv = A.inverse();
        double normAinv = Ainv.normInf();
        return normA*normAinv;
    }

    public static void test(Matrix A, mVector b, int p) throws Exception {
        System.out.println("n = " + b.getLength());
        System.out.println("Число обусловленности матрицы A : " + cond(A));
        mVector x = gMRESwithILUcond(A,b,p);
        System.out.println("_________Vector X__________");
        //x.print();
        Matrix res = A.solve(b.toMatrix());
        for (int i = 0; i < x.getLength();i++){
            System.out.print(x.get(i));
            System.out.print("  =  ");
            System.out.print(res.get(i,0));
            System.out.println();
        }
        System.out.println("=========================================================");

    }
    public static void main(String[] args) throws Exception {
        int p = 3;
        System.out.println("ХОРОШО ОБУСЛОВЛЕННЫЕ МАТРИЦЫ!");

        //---------------------------test 1-----------------------------------------//
        int n = 10;
        double [][] matrixA = {
                { 4372.0, -2.0, -4.0, 0.0, -5.0, 0.0, -3.0, -4.0, -1.0, -3.0},
                { -4.0, 3011.0, -3.0, -3.0, -4.0, -3.0, 0.0, -2.0, -1.0, -2.0},
                { -1.0, -3.0, 4114.0, -4.0, -2.0, -4.0, -3.0, -1.0, -3.0, -1.0},
                {-5.0, 0.0, -3.0, 4686.0, 0.0, -1.0, -2.0, 0.0, -2.0, -2.0},
                {0.0, -1.0, -5.0, -5.0, 4656.0, -4.0, -4.0, 0.0, -5.0, -5.0},
                { -2.0, 0.0, -2.0, -4.0, -3.0, 4469.0, -2.0, -4.0, 0.0, -2.0},
                {0.0, -4.0, -5.0, -3.0, -4.0, 0.0, 3163.0, -4.0, -1.0, -1.0},
                {-5.0, -2.0, -5.0, -4.0, 0.0, -4.0, -2.0, 4547.0, -2.0, -3.0},
                {-1.0, -1.0, -1.0, -4.0, -4.0, -1.0, -3.0, -4.0, 3969.0, 0.0},
                {-5.0, -2.0, -4.0, 0.0, -4.0, -2.0, -3.0, -2.0, -2.0, 1162.0},
        };

        mVector b = new mVector(n);
        b.vector[0] = 3468.0;
        b.vector[1] = 221.0;
        b.vector[2] = 4422.0;
        b.vector[3] = 1165.0;
        b.vector[4] = 2183.0;
        b.vector[5] = 4504.0;
        b.vector[6] = 431.0;
        b.vector[7] = 261.0;
        b.vector[8] = 4132.0;
        b.vector[9] = 2417.0;

        /*test(new Matrix(matrixA),b,p);

        matrixA[0][0] -= 50;
        matrixA[0][4] +=5;
        matrixA[9][0] +=3;
        b.set(9,b.get(9)-100);*/

        Matrix m_matrixA = new Matrix(matrixA);
        test(m_matrixA,b,p);



        //-----------------------------test 2-------------------------------------
        n = 20;
        m_matrixA = generateMatrix(n);
        b = generateVector(n);
        test(m_matrixA,b,p);

        //-----------------------------test 3-------------------------------------
        n = 30;
        m_matrixA = generateMatrix(n);
        b = generateVector(n);
        test(m_matrixA,b,p);

        System.out.println("ПЛОХО ОБУСЛОВЛЕННЫЕ МАТРИЦЫ!");

        //-----------------------------test 4-------------------------------------
        n = 11;
        m_matrixA = generateMatrixBadCond(n);
        //printMatrix(m_matrixA);
        b = generateVector(n);
        //b.print();
        test(m_matrixA,b,p);

        //---------------------------test 5-----------------------------------------//
        n = 20;
        double [][] matrixA_ = {
                { 26.0, -4.0,    0.0, -3.0, 0.0, -5.0, 0.0, -2.0, 0.0, -3.0, 0.0, -4.0, 0.0, 0.0, 0.0, -2.0, 0.0, -5.0, 0.0, -4.0},
                { 0.0,   4668.0,  0.0, -2.0, 0.0, -5.0, 0.0, -2.0, 0.0, -1.0, 0.0, -3.0, 0.0, -3.0, 0.0, -2.0, 0.0, -3.0, 0.0, -1.0},
                { 0.0,  -1.0, 93.0, -1.0, 0.0, -1.0, 0.0, -4.0, 0.0, -2.0, 0.0, -4.0, 0.0, -1.0, 0.0, -3.0, 0.0, -1.0, 0.0, 0.0},
                { 0.0,  -3.0, 0.0, 29649.0, 0.0, -2.0, 0.0, -5.0, 0.0, -2.0, 0.0, -4.0, 0.0, 0.0, 0.0, -5.0, 0.0, -3.0, 0.0, -2.0},
                { 0.0,   0.0, 0.0, -4.0, 94.0, -2.0, 0.0, -2.0, 0.0, -4.0, 0.0, -2.0, 0.0, -4.0, 0.0, -4.0, 0.0, -2.0, 0.0, -5.0},
                { 0.0,  -4.0, 0.0, -3.0, 0.0, 11914.0, 0.0, -5.0 ,0.0, -2.0, 0.0, -2.0, 0.0, -1.0, 0.0, 0.0, 0.0, -1.0, 0.0, -1.0},
                { 0.0,  -1.0, 0.0, -5.0, 0.0, -4.0, 71.0, -4.0, 0.0, 0.0, 0.0, -4.0, 0.0, 0.0, 0.0, -3.0, 0.0, -2.0, 0.0, -1.0},
                { 0.0,  -5.0, 0.0, -4.0, 0.0, -5.0, 0.0, 18748.0, 0.0, 0.0, 0.0, -2.0, 0.0, -2.0, 0.0, -3.0, 0.0, 0.0, 0.0, -3.0},
                { 0.0,  -3.0, 0.0, -4.0, 0.0, 0.0, 0.0, -5.0, 94.0, 0.0, 0.0, -5.0, 0.0, -4.0, 0.0, -2.0, 0.0, -3.0, 0.0, -3.0},
                { 0.0,  -2.0, 0.0, -3.0, 0.0, -3.0, 0.0, -5.0, 0.0, 13164.0, 0.0, -1.0, 0.0, 0.0, 0.0 ,-2.0, 0.0, -3.0, 0.0, -2.0},
                { 0.0,  -2.0, 0.0, -1.0, 0.0, -3.0, 0.0, -2.0 ,0.0, -2.0, 76.0, 0.0, 0.0, -3.0, 0.0, -1.0, 0.0, -4.0, 0.0, -1.0},
                { 0.0,   0.0, 0.0, -2.0, 0.0, -5.0, 0.0, 0.0, 0.0, 0.0, 0.0, 7625.0, 0.0, -1.0, 0.0, -5.0, 0.0, -4.0, 0.0, -1.0},
                { 0.0,  -5.0, 0.0, -5.0, 0.0, -4.0, 0.0, 0.0, 0.0, 0.0, 0.0, -4.0 ,41.0 ,-2.0, 0.0, -4.0, 0.0, -1.0, 0.0, -1.0},
                { 0.0,  -5.0, 0.0, -4.0, 0.0, -5.0, 0.0, -4.0, 0.0, -3.0, 0.0, -1.0, 0.0, 24681.0, 0.0, -3.0, 0.0, -4.0, 0.0, 0.0},
                { 0.0,  -3.0, 0.0, -5.0, 0.0, -1.0, 0.0, -3.0, 0.0, -5.0, 0.0, -4.0, 0.0, -3.0, 1.0, -1.0, 0.0, -2.0, 0.0, -1.0},
                { 0.0,  -5.0 ,0.0, -2.0, 0.0, -5.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -2.0, 0.0, 33529.0 ,0.0, 0.0, 0.0, -2.0},
                { 0.0,  -4.0, 0.0, -3.0, 0.0, -4.0, 0.0, -2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -4.0, 11.0, -5.0, 0.0, -2.0},
                { 0.0,  -4.0, 0.0, -2.0, 0.0, -5.0, 0.0, -2.0, 0.0, -3.0, 0.0, -4.0, 0.0, 0.0, 0.0, -4.0, 0.0, 32967.0, 0.0, -1.0},
                { 0.0,  -4.0, 0.0, -2.0, 0.0, 0.0, 0.0, -4.0, 0.0, -1.0, 0.0, -4.0, 0.0, -3.0, 0.0, -3.0, 0.0, -2.0, 44.0, -4.0},
                { 0.0,   0.0, 0.0, -2.0, 0.0, -4.0, 0.0, 0.0, 0.0, -1.0, 0.0, -5.0, 0.0, -5.0, 0.0, 0.0, 0.0, -2.0, 0.0, 26599.0},

        };
        double[] b_ = {3122.0, 1854.0, 4798.0, 871.0, 2457.0, 867.0, 2374.0, 3465.0, 2732.0, 833.0, 1935.0, 4414.0, 4343.0, 1933.0, 3497.0, 3025.0, 3131.0, 3365.0, 537.0, 2104.0 };
        b =  new mVector(b_);
        /*test(new Matrix(matrixA_),b,p);
        matrixA_[3][3] -= 10;
        matrixA_[1][5] +=5;
        matrixA_[12][1] +=3;
        b.set(17,b.get(17)-10);*/
        m_matrixA = new Matrix(matrixA_);
        test(m_matrixA,b,p);

        //-----------------------------test 6-------------------------------------
        n = 31;
        m_matrixA = generateMatrixBadCond(n);
        //printMatrix(m_matrixA);
        b = generateVector(n);
        //b.print();
        test(m_matrixA,b,p);

    }
}

