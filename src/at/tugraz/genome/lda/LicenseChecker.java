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

package at.tugraz.genome.lda;

import java.io.File;
import java.util.Vector;

import javax.swing.ImageIcon;

import jlk.LicenseHandler;
import jlk.MatrixLipidDataAnalyzer2;

/**
 * 
 * @author Juergen Hartler
 * @author Leonida M. Lamp
 *
 */
public class LicenseChecker
{
	public void checkLicense() 
	{
		MatrixLipidDataAnalyzer2 licenseMatrix = new MatrixLipidDataAnalyzer2();
	  Vector<String> forbiddenkeys = new Vector<String>();
	  LicenseHandler.setForbiddenKeyMD5Hashes(forbiddenkeys);
	  LicenseHandler.setLicenseArray(licenseMatrix);
	  LicenseHandler.setApplicationIcon(new ImageIcon(getClass().getResource("/images/Delete.gif")));
	  LicenseHandler.setAppName(licenseMatrix.getApplicationName());
	  LicenseHandler.setModuleNames(new String[]{"0"});
	  LicenseHandler.setDemoModeEnabled(false);
	  LicenseHandler.setUseRoamingProfileEnabled(false);
	  File licenseFolder = new File(Settings.getLicensePath());
	  File ldaHomeFolder = new File (Settings.getLdaUserHomePath());
	  if (!ldaHomeFolder.exists()) {
	    if (!ldaHomeFolder.mkdirs()) {
	      System.err.println("Could not create mandatory directory "+ldaHomeFolder.getAbsolutePath());
	      System.exit(-1);
	    }
	    try {
	      Runtime.getRuntime().exec("attrib +h \""+ldaHomeFolder.getAbsolutePath()+"\"");
	    } catch (Exception ex) {
	      System.err.println("Could not hide directory "+ldaHomeFolder.getAbsolutePath());
	    }
	  }
	  ////this is for the demo license for the reviewers
	  if (!licenseFolder.exists()) {               
	  	//TODO: this license for reviewing is valid until the 30th of November 2025!
	  	LicenseHandler.writeLicenseStringAndUserName(licenseFolder,"mx44-sGS]-[WXu-Ob6S-^7kF-4Xnv-X65r-K]^z","Test 2.11.0");
	  }
	  LicenseHandler.checkLicenseInFolder(licenseFolder);
	}
	
	public static void showLicenseDialog()
	{
		LicenseHandler.showLicenseDialog();
	}
	
	/**
	 * Set to false when pushing to public git
	 * @return
	 */
	public static boolean isCheckLicense()
	{
		return true;
	}
	
}
