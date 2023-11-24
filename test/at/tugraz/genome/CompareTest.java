package at.tugraz.genome;

import java.util.Collections;
import java.util.Vector;

import at.tugraz.genome.lda.msn.hydroxy.parser.HydroxyEncoding;
import at.tugraz.genome.lda.msn.vos.FattyAcidVO;
import at.tugraz.genome.lda.utils.StaticUtils;

public class CompareTest
{
	public static void main(String[] args)
  {
		try
		{
//			Vector<FattyAcidVO> vos = StaticUtils.decodeFAsFromHumanReadableName("18:1(n-7)_O-18:1", new HydroxyEncoding("hydroxylationEncoding.txt"),
//					new HydroxyEncoding("lcbHydroxylationEncoding.txt"), false, null);
			Vector<FattyAcidVO> vos = StaticUtils.decodeFAsFromHumanReadableName("18:0_O-18:1(n-7)", new HydroxyEncoding("hydroxylationEncoding.txt"),
					new HydroxyEncoding("lcbHydroxylationEncoding.txt"), false, null);
			
			System.out.println("before");
			for (FattyAcidVO vo : vos)
			{
				System.out.println(vo.toString());
			}
			Collections.sort(vos);
			System.out.println("after");
			for (FattyAcidVO vo : vos)
			{
				System.out.println(vo.toString());
			}
		}
		catch (Exception ex)
		{
			System.out.println("hi");
		}
  }
	
	
}
