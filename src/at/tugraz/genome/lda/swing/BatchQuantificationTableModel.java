/* 
 * This file is part of Lipid Data Analyzer
 * Lipid Data Analyzer - Automated annotation of lipid species and their molecular structures in high-throughput data from tandem mass spectrometry
 * Copyright (c) 2017 Juergen Hartler, Andreas Ziegl, Gerhard G. Thallinger 
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

package at.tugraz.genome.lda.swing;

import java.io.File;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Vector;

import javax.swing.Icon;
import javax.swing.ImageIcon;
import javax.swing.table.AbstractTableModel;

import at.tugraz.genome.lda.vos.RawQuantificationPairVO;

/**
 * 
 * @author Juergen Hartler
 *
 */
public class BatchQuantificationTableModel extends AbstractTableModel
{

  /**
   * 
   */
  private static final long serialVersionUID = -6221307059304744458L;

  private static final Integer STATUS_ADDED = new Integer(0);

  private static final Integer STATUS_UPLOADED = new Integer(1);

  private static final Integer STATUS_ERROR = new Integer(2);

  private List<RawQuantificationPairVO> fileList;

  private Map<RawQuantificationPairVO,Integer> fileStatusMap;

  private int numFilesLeft;

  private ImageIcon uploadedIcon = new ImageIcon(getClass().getResource(
      "/images/uploaded.gif"));

  private ImageIcon errorIcon = new ImageIcon(getClass().getResource(
      "/images/error.gif"));

  public BatchQuantificationTableModel()
  {
    fileList = new ArrayList<RawQuantificationPairVO>();
    fileStatusMap = new HashMap<RawQuantificationPairVO,Integer>();
    numFilesLeft = 0;
  }

  public List<RawQuantificationPairVO> getFilesLeftToQuantify()
  {
    List<RawQuantificationPairVO> list = new ArrayList<RawQuantificationPairVO>(numFilesLeft);
    Iterator<RawQuantificationPairVO> it = fileList.iterator();
    while (it.hasNext()) {
      final RawQuantificationPairVO file = it.next();
      final Integer status = (Integer) fileStatusMap.get(file);
      if (status == STATUS_ADDED) {
        list.add(file);
      }
    }
    return list;
  }

  public int getNumFiles()
  {
    return fileList.size();
  }

  public int getNumFilesLeftToUpload()
  {
    return numFilesLeft;
  }

  public int getFileRow(RawQuantificationPairVO file)
  {
    int row = 0;
    Iterator<RawQuantificationPairVO> it = fileList.iterator();
    while (it.hasNext()) {
      if (file.equals((RawQuantificationPairVO) it.next())) {
        return row;
      }
      row++;
    }
    return -1;
  }

  public void addFiles(Vector<RawQuantificationPairVO> files)
  {
    for (RawQuantificationPairVO pairVO : files) {
//      if (!fileStatusMap.containsKey(file)) {
        fileList.add(pairVO);
        fileStatusMap.put(pairVO, STATUS_ADDED);
        numFilesLeft++;
//      }
    }
    fireTableDataChanged();
  }

/*  public void removeRows(int[] rows)
  {
    for (int i = rows.length - 1; i >= 0; i--) {
      final File file = (File) fileList.remove(rows[i]);
      final Integer status = (Integer) fileStatusMap.remove(file);
      if (status != null) {
        final long fsize = file.length();
        totalSize -= fsize;
        if (status == STATUS_ADDED) {
          totalSizeLeft -= fsize;
          numFilesLeft--;
        }
      }
    }
    fireTableDataChanged();
  }*/

  public void clearFiles()
  {
    fileList.clear();
    fireTableDataChanged();
  }

  public int getColumnCount()
  {
    return 4;
  }

  public String getColumnName(int column)
  {
    switch (column) {
      case 0:
        return " ";
      case 1:
        return "Data";
      case 2:
        return "Quant-File";
      case 3:
        return "Progress";
      default:
        return null;
    }
  }

  public int getRowCount()
  {
    return fileList != null ? fileList.size() : 0;
  }

  public Object getValueAt(int row, int column)
  {
    switch (column) {
      case 0:
        final RawQuantificationPairVO file = (RawQuantificationPairVO) fileList.get(row);
        final Integer status = (Integer) fileStatusMap.get(file);
        if (status == STATUS_UPLOADED) {
          return uploadedIcon;
        } else if (status == STATUS_ERROR) {
          return errorIcon;
        } else {
          return null;
        }
      case 1:
        return ((RawQuantificationPairVO) fileList.get(row)).getRawFileName();
      case 2:
        return ((RawQuantificationPairVO) fileList.get(row)).getQuantFileName();
      case 3:
        return ((RawQuantificationPairVO) fileList.get(row)).getStatus();
      default:
        return null;
    }
  }

  @SuppressWarnings({ "unchecked", "rawtypes" })
  public Class getColumnClass(int column)
  {
    switch (column) {
      case 0:
        return Icon.class;
      case 1:
        return String.class;
      case 2:
        return String.class;
      case 3:
        return String.class;
      default:
        return null;
    }
  }

  public boolean isCellEditable(int row, int column)
  {
    return false;
  }

  public void fileQuantified(RawQuantificationPairVO file)
  {
    final Integer status = (Integer) fileStatusMap.get(file);
    if (status != STATUS_UPLOADED) {
      fileStatusMap.put(file, STATUS_UPLOADED);
      fireTableDataChanged();
      numFilesLeft--;
    }
  }

  public void fileQuantificationError(RawQuantificationPairVO file)
  {
    final Integer status = (Integer) fileStatusMap.get(file);
    if (status != STATUS_ERROR) {
      fileStatusMap.put(file, STATUS_ERROR);
      fireTableDataChanged();
      numFilesLeft--;
    }
  }
  
  /**
   * updates the raw file part of the table entry
   * @param file the quantification pair in its old form
   * @param newRawFile the new raw file
   */
  public void updateRawFile(RawQuantificationPairVO file, File newRawFile){
    final Integer status = (Integer) fileStatusMap.get(file);
    fileStatusMap.remove(file);
    file.setRawFile(newRawFile);
    fileStatusMap.put(file, status);
  }
  
  public RawQuantificationPairVO getDataByRow(int row){
    return (RawQuantificationPairVO)this.fileList.get(row);
  }
}