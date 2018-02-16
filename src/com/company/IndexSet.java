package com.company;

import Jama.Matrix;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map;

/**
 * Created by User on 07.01.2018.
 */
public class IndexSet{
    public Map<ArrayList<Integer>,Double> indexSet;

    public IndexSet(Matrix matrix) {
        this.indexSet = new HashMap<>();
        for (int i=0; i < matrix.getRowDimension(); i++){
            for (int j=0; j < matrix.getColumnDimension(); j++){
                if (matrix.get(i,j) != 0){
                    Pair indexPair = new Pair(i,j);
                    this.indexSet.put(indexPair.getPair(), matrix.get(i,j));
                }
            }
        }
    }
    public void addPortraitMatrix(Matrix matrix){
        for (int i=0; i < matrix.getRowDimension(); i++){
            for (int j=0; j < matrix.getColumnDimension(); j++){
                if (matrix.get(i,j) != 0){
                    Pair indexPair = new Pair(i,j);
                    this.indexSet.put(indexPair.getPair(), matrix.get(i,j));
                }
            }

        }
    }
    public boolean containsPair(int i,int j){
        Pair indexPair = new Pair(i,j);
        return this.indexSet.containsKey(indexPair.getPair());
    }

}