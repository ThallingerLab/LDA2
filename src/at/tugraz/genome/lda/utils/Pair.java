/* 
 * This file is part of Lipid Data Analyzer
 * Lipid Data Analyzer - Automated annotation of lipid species and their molecular structures in high-throughput data from tandem mass spectrometry
 * Copyright (c) 2017 Juergen Hartler, Andreas Ziegl, Gerhard G. Thallinger, Leonida M. Lamp 
 * DO NOT ALTER OR REMOVE COPYRIGHT NOTICES OR THIS FILE HEADER. 
 *  
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * by the Free Software Foundation, either version 3 of the License, or 
 * (at your option) any later version.
 *  
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details. 
 *  
 * You should have received a copy of the GNU General Public License 
 * along with this program. If not, see <http://www.gnu.org/licenses/>.
 *
 * Please contact lda@genome.tugraz.at if you need additional information or 
 * have any questions.
 */ 
package at.tugraz.genome.lda.utils;
 
/**
 * Pair Class using generics
 * 
 * @author Leonida M. Lamp
 *
 */
public class Pair<K, V>
{
    public final K key;
    public final V value;
 
    /**
     * Constructor for a Pair
     * @param key
     * @param value
     */
    public Pair(K key, V value)
    {
        this.key = key;
        this.value = value;
    }
    
    public K getKey() {
      return this.key;
    }
    
    public V getValue() {
      return this.value;
    }
 
    @Override
    public boolean equals(Object o)
    {
        if (this == o) {
            return true;
        }
        if (o == null || getClass() != o.getClass()) {
            return false;
        }
        Pair<?, ?> pair = (Pair<?, ?>) o;
        if (!key.equals(pair.key)) {
            return false;
        }
        return value.equals(pair.value);
    }
 
    /**
     * Computes hash code for an object to support hash tables
     * using hash codes of the underlying objects
     */
    @Override
    public int hashCode()
    {
        return 31 * key.hashCode() + value.hashCode();
    }
 
    @Override
    public String toString() {
        return "(" + key + ", " + value + ")";
    }
    
    /**
     * Factory method for creating an immutable instance of a typed Pair
     */
    public static <K, V> Pair <K, V> of(K a, V b)
    {
        return new Pair<>(a, b);
    }
}
