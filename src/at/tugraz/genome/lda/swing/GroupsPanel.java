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

import java.awt.Dimension;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.Insets;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.File;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;
import java.util.Hashtable;
import java.util.List;
import java.util.Set;
import java.util.Vector;

import javax.swing.JFrame;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JTable;

import at.tugraz.genome.lda.WarningMessage;
import at.tugraz.genome.lda.utils.StaticUtils;

/**
 * 
 * @author Juergen Hartler
 *
 */
public class GroupsPanel extends JPanel implements ActionListener
{

  private static final long serialVersionUID = 208895624290731494L;
  
  private Vector<String> groups_ = new Vector<String>();
  private Hashtable<String,Vector<File>> filesOfGroup_ = new Hashtable<String,Vector<File>>();
  private Hashtable<String,GroupPanel> panelOfGroup_ = new Hashtable<String,GroupPanel>();
  private final int itemsInRow_ = 3;
  private int maxWidth_;
  private final int scrollBarWidth_ = 10;
 
  public GroupsPanel(){
    super();
    setLayout(new GridBagLayout());
  }
  
  public void addGroup(String groupName, Vector<File> selectedFiles, int maxWidth){
    groups_.add(groupName);
    filesOfGroup_.put(groupName, selectedFiles);
    maxWidth_ = maxWidth;
    this.visualizeGroups();
  }
  
  private void visualizeGroups(){
    this.removeAll();
    int count = 0;
    int insetLeft = 5;
    int insetRight = 5;
    JPanel panelToAddGroupPanel = this;
    if (groups_.size()>itemsInRow_){
      panelToAddGroupPanel = new JPanel();
    }
    for (String group : groups_){
      GroupPanel groupPanel = new GroupPanel(group,filesOfGroup_.get(group),maxWidth_/itemsInRow_-insetLeft-insetRight-scrollBarWidth_,this);
      panelOfGroup_.put(group, groupPanel);
      panelToAddGroupPanel.add(groupPanel,new GridBagConstraints(count%itemsInRow_, count/itemsInRow_, 1, 1, 0.0, 0.0
          ,GridBagConstraints.CENTER, GridBagConstraints.NONE, new Insets(2, insetLeft, 2, insetRight), 0, 0));
      count++;
    }
    if (groups_.size()>itemsInRow_){
      panelToAddGroupPanel.setPreferredSize(new Dimension(maxWidth_-insetLeft-insetRight,(count/itemsInRow_+1)*150));
      JScrollPane scrollPane = new JScrollPane(panelToAddGroupPanel);
      scrollPane.setPreferredSize(new Dimension(maxWidth_-insetLeft-insetRight-scrollBarWidth_, 170));
      this.add(scrollPane);
    }

    this.invalidate();
    this.updateUI();
  }
  
  public Vector<String> getGroups(){
    return this.groups_;
  }
  
  public void removeAllGroups(){
    this.removeAll();
    groups_ = new Vector<String>();
    filesOfGroup_ = new Hashtable<String,Vector<File>>();
    panelOfGroup_ = new Hashtable<String,GroupPanel>();
    this.invalidate();
    this.updateUI();
  }
  
  public void removeFiles(Vector<File> filesToRemove){
    Vector<String> newGroups = new Vector<String>();
    for (String group:groups_){
      Vector<File> newFiles = new Vector<File>();
      for (File file:filesOfGroup_.get(group)){
        boolean remove = false;
        for (File removeFile:filesToRemove){
          if (file.getAbsolutePath().equalsIgnoreCase(removeFile.getAbsolutePath()))
            remove = true;
        }
        if (!remove)
          newFiles.add(file);
      }
      if (newFiles.size()>0){
        filesOfGroup_.put(group,newFiles);
        newGroups.add(group);
      }else{
        filesOfGroup_.remove(group);
        panelOfGroup_.remove(group);
      }  
    }
    this.groups_ = newGroups;
    this.visualizeGroups();
  }
  
  public void actionPerformed(ActionEvent e)
  {
    String command = e.getActionCommand();
    if (command.startsWith("renameGroup;")){
      String groupName = command.substring("renameGroup;".length());
      InputDialog dlg = new InputDialog(new JFrame(), "Rename group", "Enter the group name", String.valueOf(groupName));
      String newGroupName = dlg.getEnteredText();
      if (newGroupName!=null&&newGroupName.length()>0){
        newGroupName = newGroupName.trim();
        Vector<File> selectedFiles = filesOfGroup_.get(groupName);
        filesOfGroup_.remove(groupName);
        filesOfGroup_.put(newGroupName, selectedFiles);
        panelOfGroup_.put(newGroupName, panelOfGroup_.get(groupName));
        panelOfGroup_.remove(groupName);
        Vector<String> newGroups = new Vector<String>();
        for (String group : groups_){
          if (group.equalsIgnoreCase(groupName))
            newGroups.add(newGroupName);
          else
            newGroups.add(group);
        }
        groups_ = newGroups;
        this.visualizeGroups();
      }else{
        new WarningMessage(new JFrame(), "Warning", "You have to specifiy a name for the group!");
      }
    }
    if (command.startsWith("removeItems;")){
      String groupName = command.substring("removeItems;".length());
      GroupPanel panel = this.panelOfGroup_.get(groupName);
      JTable table = panel.getTable();
      int[] selectedColumns = table.getSelectedRows();
      if (selectedColumns!=null && selectedColumns.length>0){
        List<Integer> selectedList = new ArrayList<Integer>();
        for (int i=0;i!=selectedColumns.length;i++){
          selectedList.add(selectedColumns[i]);
        }
        Collections.sort(selectedList);
//        Vector<File> filesToRemove = new Vector<File>();
        Vector<File> filesOfGroup = filesOfGroup_.get(groupName);
        for (int i=(selectedList.size()-1);i!=-1;i--){
          filesOfGroup.remove((int)selectedList.get(i));
        }
        if (filesOfGroup.size()>0){
          filesOfGroup_.put(groupName,filesOfGroup);
          panel.updateFiles(filesOfGroup);
        }else{
          command = "deleteGroup;"+groupName;
        }
      }
    }
    if (command.startsWith("deleteGroup;")){
      String groupName = command.substring("deleteGroup;".length());
      for (int i=0; i!=groups_.size();i++){
        if (groups_.get(i).equalsIgnoreCase(groupName)){
          filesOfGroup_.remove(groupName);
          panelOfGroup_.remove(groupName);
          groups_.remove(i);
          break;
        }
      }
      this.visualizeGroups();
    }    
  }
  
  public Set<String> getExpsOfGroupOneExpBelongsTo(String expirementName, Vector<String> allExperimentNames){
    Vector<String> groupsTheExpBelongsTo = new Vector<String>();
    Hashtable<String,Vector<File>> groupFiles = filesOfGroup_;
    for (String group : groups_){
      boolean belongsToGroup = false;
      for (File file : groupFiles.get(group)){
        String fileName = StaticUtils.extractFileName(file.getAbsolutePath());
        if (fileName.contains(expirementName))
          belongsToGroup = true;
      }
      if (belongsToGroup)
        groupsTheExpBelongsTo.add(group);
    }
    if (groupsTheExpBelongsTo.size()>0){
      Hashtable<String,String> expsUpdate = new Hashtable<String,String>();
      for (String group:groupsTheExpBelongsTo){
        for (File file : groupFiles.get(group)){
          String fileName = StaticUtils.extractFileName(file.getAbsolutePath());
          for (String expName : allExperimentNames){
            if (fileName.contains(expName)){
              expsUpdate.put(expName, expName);
              break;
            }  
          }
        }
      }
      return expsUpdate.keySet();
    }
    return new HashSet<String>();
  }
  
  public Hashtable<String,Vector<File>> getGroupFiles(){
    return this.filesOfGroup_;
  }
}
