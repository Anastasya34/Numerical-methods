package com.company;

import Jama.Matrix;

/**
 * Created by User on 29.01.2018.
 */
public class ILUResult {
    Matrix L;
    Matrix U;

    public ILUResult(Matrix L, Matrix U) {
        this.L = L;
        this.U = U;
    }

    public Matrix getL() {
        return L;
    }
    public void setL(Matrix l) {
        L = l;
    }


    public void setU(Matrix u) {
        U = u;
    }

    public Matrix getU() {
        return U;
    }
    public Matrix getLU() {
        return L.times(U);
    }

}
