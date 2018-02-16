package com.company;

import java.util.ArrayList;

/**
 * Created by User on 07.01.2018.
 */
public class Pair{
    int i,j;

    public Pair(int i, int j) {
        this.i = i;
        this.j = j;
    }

    public ArrayList<Integer> getPair(){
        ArrayList<Integer> indexPair = new ArrayList<>();
        indexPair.add(i);
        indexPair.add(j);
        return indexPair;
    }
}