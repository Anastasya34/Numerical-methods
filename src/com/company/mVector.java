package com.company;

import Jama.Matrix;

/**
 * Created by User on 07.01.2018.
 */
public class mVector {
    public double[] vector;
    private int length;

    public mVector(int length) throws Exception {
        if (length < 1) {
            throw new Exception("Неверный размер.");
        }
        this.length = length;
        vector = new double[length];
    }
    public mVector(int length, double elem) throws Exception {
        if (length < 1) {
            throw new Exception("Неверный размер.");
        }
        this.length = length;
        vector = new double[length];
        for (int i = 0; i < length; i++){
            vector[i] = elem;
        }
    }
    public mVector(double[] v) throws Exception {
        this.length = v.length;
        this.vector = v;
    }
    public mVector(Matrix A, int column) throws Exception {
        this.length = A.getRowDimension();
        vector = new double[length];
        for (int j = 0; j < length; j++) {
            this.vector[j] = A.get(j,column);
        }
    }

    public int getLength() {
        return length;
    }

    public void print() {
        for (double item : vector) {
            System.out.println(item + " ");
        }
    }
    public mVector getVector(int start,int stop) throws Exception {
        if (stop >= this.getLength() || start >= this.getLength() ) {
            throw new Exception("Выход за пределы размерности вектора.");
        }
        mVector v = new mVector(stop-start+1);
        for (int j=0, i = start; i <= stop; i++, j++){
            v.set(j, this.get(i));
        }
        return v;
    }

    public double mul(mVector second) throws Exception {
        if (length != second.getLength()) {
            throw new Exception("Неверный вектор.");
        }
        double result = 0;
        for (int i = 0; i < length; i++) {
            result += this.vector[i] * second.vector[i];
        }
        return result;
    }


    public mVector mul(double num) throws Exception {
        mVector result = new mVector(length);
        for (int i = 0; i < length; i++) {
            result.vector[i] = this.vector[i] * num;
        }
        return result;
    }
    public mVector div(double num) throws Exception {
        mVector result = new mVector(length);
        for (int i = 0; i < length; i++) {
            result.vector[i] = this.vector[i] / num;
        }
        return result;
    }
    public mVector add(mVector sub) throws Exception {
        if (length != sub.getLength()) {
            throw new Exception("Неверный вектор.");
        }
        mVector result = new mVector(length);
        for (int i = 0; i < length; i++) {
            result.vector[i] = this.vector[i] + sub.vector[i];
        }
        return result;
    }
    public mVector subtract(mVector sub) throws Exception {
        if (length != sub.getLength()) {
            throw new Exception("Неверный вектор.");
        }
        mVector result = new mVector(length);
        for (int i = 0; i < length; i++) {
            result.vector[i] = this.vector[i] - sub.vector[i];
        }
        return result;
    }
    public double norm(){
        double a = 0.0;
        for(int i = 0; i < vector.length; i++){
            a += vector[i]*vector[i];
        }
        return Math.sqrt(a);
    }

    public double normI() {
        double max = Math.abs(vector[0]);
        for (int i = 1; i < length; i++) {
            if (Math.abs(vector[i]) > max) {
                max = Math.abs(vector[i]);
            }
        }
        return max;
    }
    public Matrix toMatrix(){
        double [][] matrix = new double[this.getLength()][1];
        for (int i = 0; i < this.getLength(); i++){
            matrix [i][0] = this.get(i);
        }
        return new Matrix(matrix);
    }

    public double get(int i) {
        return vector[i];
    }

    public void set(int i,double elem) {
        this.vector[i] = elem;
    }
}
