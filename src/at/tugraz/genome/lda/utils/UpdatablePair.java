/* 
 * This file is part of Lipid Data Analyzer
 * Lipid Data Analyzer - Automated annotation of lipid species and their molecular structures in high-throughput data from tandem mass spectrometry
 * Copyright (c) 2023 Juergen Hartler, Andreas Ziegl, Gerhard G. Thallinger, Leonida M. Lamp 
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
 * Pair class that is updatable in contrast to the otherwise analogous Pair class offered by javafx.
 * 
 * @author Leonida M. Lamp
 *
 */
public class UpdatablePair<K, V>
{
  private K key_;
  private V value_;
 
  /**
   * Constructor for a Pair
   * @param key
   * @param value
   */
  public UpdatablePair(K key, V value)
  {
    this.key_ = key;
    this.value_ = value;
  }
    
  public K getKey() 
  {
    return this.key_;
  }
    
  public V getValue() 
  {
    return this.value_;
  }
    
  public void setKey(K key)
  {
  	this.key_ = key;
  }
    
  public void setValue(V value)
  {
  	this.value_ = value;
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
    UpdatablePair<?, ?> pair = (UpdatablePair<?, ?>) o;
    if (!key_.equals(pair.key_)) {
      return false;
    }
    return value_.equals(pair.value_);
  }
 
  /**
   * Computes hash code for an object to support hash tables
   * using hash codes of the underlying objects
   */
  @Override
  public int hashCode()
  {
    return 31 * key_.hashCode() + value_.hashCode();
  }
 
  @Override
  public String toString() {
    return "(" + key_ + ", " + value_ + ")";
  }
    
  /**
   * Factory method for creating an immutable instance of a typed Pair
   */
  public static <K, V> UpdatablePair <K, V> of(K a, V b)
  {
    return new UpdatablePair<>(a, b);
  }
}
