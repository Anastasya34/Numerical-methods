package com.company;

import Jama.Matrix;

/**
 * Created by User on 25.01.2018.
 */
public class ArnoldiResult {
    Matrix H;
    Matrix V;

    public ArnoldiResult(Matrix h, Matrix v) {
        H = h;
        V = v;
    }

    public Matrix getH() {
        return H;
    }

    public void setH(Matrix h) {
        H = h;
    }

    public Matrix getV() {
        return V;
    }

    public void setV(Matrix v) {
        V = v;
    }
}
