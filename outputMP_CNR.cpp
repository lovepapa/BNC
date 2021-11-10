//<<test BDS MP and CNR  2017-011-15
//<加入跳变探测
#pragma once
#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <iomanip>
#include "newmat.h"
#include "rnxobsfile.h"
#include "rnxnavfile.h"
#include "ephemeris.h"//<星历数据结构
#include "public.h"//<观测量数据结构
#include "bancroft.h"
#include "panutils.h"

using namespace std;

int main(int argc, char *argv[])
{
	struct filterData
	{
			
		filterData(): wsize(1), mean(.0){};
				
		int wsize;
		
		double mean;
	};
	std::map<int, filterData>B1I_data, B1C_data, B1A_data, B2a_data, B2b_data, B3I_data, B3Q_data, B3A_data, LMP_data;
	//<           <arc  data>
	//----------------------------------------------------------------------------
	//<读入配置文件参数
	string navfile, obsfile, ranfile, obstype;
	char satsys('C');
	//int prn(0);
	ColumnVector TruePos(3);

	fstream cfgfile;
	cfgfile.open(".\\txy0112\\bnctxy2.conf",ios::in);
	string line;
	while ( getline(cfgfile, line) ){		 
		string key(line, 0, 30);
		strip(key);
		if (key =="obs file"){
			//obsfile = line.substr(30,20);
			obsfile = line.substr(30, 30);
		}
		else if(key =="nav file"){
			//navfile = line.substr(30,20);
			navfile = line.substr(30, 30);
		}
		//else if(key =="ran file"){
		//	ranfile = line.substr(30,20);
		//}
		else if(key =="system"){
			string temp =  line.substr(30,20);
			satsys = temp[0];
		} 
		//else if(key =="prn"){         // conf 文件指定输出一颗卫星的变量
		//	prn = asInt(line.substr(30,20));
		//}
		//else if(key =="obstype"){
		//	obstype = line.substr(30,20);
		//}
		else if(key =="true pos x"){
			TruePos(1)  = asDouble( line.substr(30,20) );
		}
		else if(key =="true pos y"){
			TruePos(2)  = asDouble( line.substr(30,20) );
		}
		else if(key =="true pos z"){
			TruePos(3) = asDouble( line.substr(30,20) );
		}
		else if(key =="end of conf"){
			break;
		}

	}

		//----------------------------------------------------------------------------
		//<输出文件
		fstream outfile(".\\txy0112\\newresult_2\\TJH21200.21mp", ios::out);
		//fstream outfile1(".\\txy0112\\check.out", ios::out);
		outfile.precision(3); outfile.setf(ios::fixed);
		/* outfile << setw(23) << "" << "SAT" << setw(11) << "AZI" << setw(9) << "ELE"
			//<< setiosflags(ios::left) << setw(8) << "mpB1I" << setiosflags(ios::left) << setw(8) << "mpB1C" << setiosflags(ios::left) << setw(8) << "mpB1A" << setiosflags(ios::left) << setw(8) << "mpB2a" << setiosflags(ios::left) << setw(8) << "mpB2b" << setiosflags(ios::left) << setw(8) << "mpB3I" << setiosflags(ios::left) << setw(8) << "mpB3Q" << setiosflags(ios::left) << setw(8) << "mpB3A"
			//<< setiosflags(ios::left) << setw(14) << "S1I" << setiosflags(ios::left) << setw(14) << "S1C" << setiosflags(ios::left) << setw(14) << "S1A" << setiosflags(ios::left) << setw(14) << "S2a" << setiosflags(ios::left) << setw(14) << "S2b" << setiosflags(ios::left) << setw(14) << "S3I" << setiosflags(ios::left) << setw(14) << "S3Q" << setiosflags(ios::left) << setw(14) << "S3A" << endl;
			<< setiosflags(ios::right) << setw(8) << "S1I" << setiosflags(ios::right) << setw(9) << "mpB1I"
			<< setiosflags(ios::right) << setw(8) << "mpB2a" << setiosflags(ios::right) << setw(8) << "mpB2b"
			<< setiosflags(ios::right) << setw(8) << "S3I" << setiosflags(ios::right) << setw(9) << "mpB3I"
			<< setiosflags(ios::right) << setw(8) << "GFIF" << endl; 输出显示*/
	for (int prn = 1; prn <=60; prn++) 
	{
		//-----------read navigation file----------------------------------
		t_rnxNavFile rnxeph;
		rnxeph.readbody(navfile);
		t_ephBDS bdseph;

		//-----------read obs file----------------------------------	
		t_rnxObsFile grin, back_grin;
		grin.readbody(obsfile);
		back_grin.readbody(obsfile);


		//-----------end of file----------------------------------	

		//<找出观测量在数据结构中的位置
		int idx_B1I(0), idx_B1C(0), idx_B1A(0), idx_B2a(0), idx_B2b(0), idx_B2(0), idx_B3I(0), idx_B3Q(0), idx_B3A(0);
		int idx_L1I(0), idx_L1C(0), idx_L1A(0), idx_L2a(0), idx_L2b(0), idx_L2(0), idx_L3I(0), idx_L3Q(0), idx_L3A(0);
		int idx_S1I(0), idx_S1C(0), idx_S1A(0), idx_S2a(0), idx_S2b(0), idx_S2(0), idx_S3I(0), idx_S3Q(0), idx_S3A(0);
		int idx_B1X(0), idx_L1X(0), idx_S1X(0);
		//<rinex version 3.03
		for (int i = 0; i < grin.nTypes(satsys); i++) {
			//if (grin.obsType(satsys, i) == "C1I")   idx_B1I = i;
			if (grin.obsType(satsys, i) == "C1X")   idx_B1X = i; //URUM站 BDS3 B1
			if (grin.obsType(satsys, i) == "C2I")   idx_B1I = i; //URUM站 BDS2/3
			if (grin.obsType(satsys, i) == "C1D")   idx_B1C = i;
			if (grin.obsType(satsys, i) == "C1A") idx_B1A = i;

			if (grin.obsType(satsys, i) == "C5D")   idx_B2a = i; //BDS3 txy 
			if (grin.obsType(satsys, i) == "C7I") idx_B2b = i; //BDS3 txy 
			if (grin.obsType(satsys, i) == "C7")  idx_B2 = i;
			if (grin.obsType(satsys, i) == "C6I")  idx_B3I = i;
			if (grin.obsType(satsys, i) == "C6Q") idx_B3Q = i;
			if (grin.obsType(satsys, i) == "C6A") idx_B3A = i;

			
			//if (grin.obsType(satsys, i) == "L1I")   idx_L1I = i;
			if (grin.obsType(satsys, i) == "L1X")   idx_L1X = i;
			if (grin.obsType(satsys, i) == "L2I")   idx_L1I = i; //
			if (grin.obsType(satsys, i) == "L1D")   idx_L1C = i;
			if (grin.obsType(satsys, i) == "L1A") idx_L1A = i;
			if (grin.obsType(satsys, i) == "L5D")   idx_L2a = i;    //BDS3 txy 
			if (grin.obsType(satsys, i) == "L7I") idx_L2b = i;       //BDS2 txy 
			if (grin.obsType(satsys, i) == "L7")  idx_L2 = i;
			if (grin.obsType(satsys, i) == "L6I")  idx_L3I = i;
			if (grin.obsType(satsys, i) == "L6Q")  idx_L3Q = i;
			if (grin.obsType(satsys, i) == "L6A") idx_L3A = i;

			//if (grin.obsType(satsys, i) == "S1I")   idx_S1I = i;
			if (grin.obsType(satsys, i) == "S1X")   idx_S1X = i;
			if (grin.obsType(satsys, i) == "S2I")   idx_S1I = i;
			if (grin.obsType(satsys, i) == "S1D")   idx_S1C = i;
			if (grin.obsType(satsys, i) == "S1A") idx_S1A = i;
			if (grin.obsType(satsys, i) == "S5D")   idx_S2a = i;
			if (grin.obsType(satsys, i) == "S7I") idx_S2b = i;
			if (grin.obsType(satsys, i) == "S7")  idx_S2 = i;
			if (grin.obsType(satsys, i) == "S6I")  idx_S3I = i;
			if (grin.obsType(satsys, i) == "S6Q") idx_S3Q = i;
			if (grin.obsType(satsys, i) == "S6A") idx_S3A = i;
		}


		//-----------Start process----------------------------------	

		double B1I(.0), B1C(.0), B1A(.0), B2a(.0), B2b(.0), B2(.0), B3I(.0), B3Q(.0), B3A(.0);
		double L1I(.0), L1C(.0), L1A(.0), L2a(.0), L2b(.0), L2(.0), L3I(.0), L3Q(.0), L3A(.0);
		double S1I(.0), S1C(.0), S1A(.0), S2a(.0), S2b(.0), S2(.0), S3I(.0), S3Q(.0), S3A(.0);
		double MP_B1I(.0), MP_B1C(.0), MP_B1A(.0), MP_B2a(.0), MP_B2b(.0), MP_B2(.0), MP_B3I(.0), MP_B3Q(.0), MP_B3A(.0), LMP(.0);
		double pre_B1I(.0), pre_B1C(.0), pre_B1A(.0), pre_B2a(.0), pre_B2b(.0), pre_B2(.0), pre_B3I(.0), pre_B3Q(.0), pre_B3A(.0), pre_LMP(.0);
		double dL1IL2a(.0), dL1IL2b(.0), pre_dL1IL2a(.0), pre_dL1IL2b(.0);
		int B1I_arc(0), B1C_arc(0), B1A_arc(0), B2a_arc(0), B2b_arc(0), B3I_arc(0), B3Q_arc(0), B3A_arc(0), LMP_arc(0);
		int L1I_arc(0), L1C_arc(0), L1A_arc(0), L2a_arc(0), L2b_arc(0), L3I_arc(0), L3Q_arc(0), L3A_arc(0);
		bncTime curtime;
		const t_rnxEpo* rnxEpo;

		double c(t_CST::freq(t_frequency::C2, 0) / (t_CST::freq(t_frequency::C2, 0) + t_CST::freq(t_frequency::C7, 0)));//B1I B2b
		double c1(t_CST::freq(t_frequency::C2, 0) / (t_CST::freq(t_frequency::C2, 0) + t_CST::freq(t_frequency::E5, 0)));//B1I B2a B
		double c2(t_CST::freq(t_frequency::E1, 0) / (t_CST::freq(t_frequency::E1, 0) + t_CST::freq(t_frequency::C7, 0)));//B1C B2b
		double c3(t_CST::freq(t_frequency::E1, 0) / (t_CST::freq(t_frequency::E1, 0) + t_CST::freq(t_frequency::E5, 0)));//B1C B2a 

		double c4(t_CST::freq(t_frequency::C2, 0) / (t_CST::freq(t_frequency::C2, 0) + t_CST::freq(t_frequency::C6, 0)));//B1I B3I
		double c5(t_CST::freq(t_frequency::C2, 0) / (t_CST::freq(t_frequency::C2, 0) + t_CST::freq(t_frequency::E8, 0)));//B1I B2
		double c6(t_CST::freq(t_frequency::E1, 0) / (t_CST::freq(t_frequency::E1, 0) + t_CST::freq(t_frequency::C6, 0)));


		double d(t_CST::freq(t_frequency::C7, 0) / (t_CST::freq(t_frequency::C2, 0) + t_CST::freq(t_frequency::C7, 0)));//B1I B2b
		double d1(t_CST::freq(t_frequency::E5, 0) / (t_CST::freq(t_frequency::C2, 0) + t_CST::freq(t_frequency::E5, 0)));//B1I B2a
		double d2(t_CST::freq(t_frequency::C7, 0) / (t_CST::freq(t_frequency::E1, 0) + t_CST::freq(t_frequency::C7, 0)));//B1C B2b
		double d3(t_CST::freq(t_frequency::E5, 0) / (t_CST::freq(t_frequency::E1, 0) + t_CST::freq(t_frequency::E5, 0)));//B1C B2a

		double d4(t_CST::freq(t_frequency::C6, 0) / (t_CST::freq(t_frequency::C2, 0) + t_CST::freq(t_frequency::C6, 0)));//B1I B3I
		double d5(t_CST::freq(t_frequency::E8, 0) / (t_CST::freq(t_frequency::C2, 0) + t_CST::freq(t_frequency::E8, 0)));//B1I B2
		double d6(t_CST::freq(t_frequency::C7, 0) / (t_CST::freq(t_frequency::E1, 0) + t_CST::freq(t_frequency::C7, 0)));//
		double d7(t_CST::freq(t_frequency::C6, 0) / (t_CST::freq(t_frequency::E1, 0) + t_CST::freq(t_frequency::C6, 0)));//B1C B3I

		double e(t_CST::freq(t_frequency::C2, 0) / (t_CST::freq(t_frequency::C2, 0) - t_CST::freq(t_frequency::C7, 0)));//B1I B2b
		double e1(t_CST::freq(t_frequency::C2, 0) / (t_CST::freq(t_frequency::C2, 0) - t_CST::freq(t_frequency::E5, 0)));//B1I B2a
		double e2(t_CST::freq(t_frequency::E1, 0) / (t_CST::freq(t_frequency::E1, 0) - t_CST::freq(t_frequency::C7, 0)));//B1C B2b
		double e3(t_CST::freq(t_frequency::E1, 0) / (t_CST::freq(t_frequency::E1, 0) - t_CST::freq(t_frequency::E5, 0)));//B1C B2a

		double e4(t_CST::freq(t_frequency::C2, 0) / (t_CST::freq(t_frequency::C2, 0) - t_CST::freq(t_frequency::C6, 0)));//B1I B3I
		double e5(t_CST::freq(t_frequency::C2, 0) / (t_CST::freq(t_frequency::C2, 0) - t_CST::freq(t_frequency::E8, 0)));//B1I B2
		double e6(t_CST::freq(t_frequency::E1, 0) / (t_CST::freq(t_frequency::E1, 0) - t_CST::freq(t_frequency::C6, 0)));//

		double f(t_CST::freq(t_frequency::C7, 0) / (t_CST::freq(t_frequency::C2, 0) - t_CST::freq(t_frequency::C7, 0)));//B1I B2b
		double f1(t_CST::freq(t_frequency::E5, 0) / (t_CST::freq(t_frequency::C2, 0) - t_CST::freq(t_frequency::E5, 0)));//B1I B2a
		double f2(t_CST::freq(t_frequency::C7, 0) / (t_CST::freq(t_frequency::E1, 0) - t_CST::freq(t_frequency::C7, 0)));//B1C B2b
		double f3(t_CST::freq(t_frequency::E5, 0) / (t_CST::freq(t_frequency::E1, 0) - t_CST::freq(t_frequency::E5, 0)));//B1C B2a

		double f4(t_CST::freq(t_frequency::C6, 0) / (t_CST::freq(t_frequency::C2, 0) - t_CST::freq(t_frequency::C6, 0)));//B1I B3I
		double f5(t_CST::freq(t_frequency::E8, 0) / (t_CST::freq(t_frequency::C2, 0) - t_CST::freq(t_frequency::E8, 0)));//B1I B2
		double f6(t_CST::freq(t_frequency::C7, 0) / (t_CST::freq(t_frequency::E1, 0) - t_CST::freq(t_frequency::C7, 0)));//
		double f7(t_CST::freq(t_frequency::C6, 0) / (t_CST::freq(t_frequency::E1, 0) - t_CST::freq(t_frequency::C6, 0)));//B1C B3I
																										 //<process every epoch
		while (rnxEpo = grin.nextEpoch()) 
		{

			for (int i = 0; i < rnxEpo->rnxSat.size(); i++) //<遍历各颗卫星，取得相应观测量	
			{		
				const t_rnxSat& rnxSat = rnxEpo->rnxSat[i];
				/*	prn = rnxSat.satNum;*/
				if (rnxSat.satSys == satsys && rnxSat.satNum == prn) {


					B1I = rnxSat.obs[idx_B1I];
					if (B1I == 0) continue;

					B1C = rnxSat.obs[idx_B1C];
					B1A = rnxSat.obs[idx_B1A];

					B2a = rnxSat.obs[idx_B2a];
					B2b = rnxSat.obs[idx_B2b];
					B2 = rnxSat.obs[idx_B2];

					B3I = rnxSat.obs[idx_B3I];
					if (B3I == 0) continue;

					B3Q = rnxSat.obs[idx_B3Q];
					B3A = rnxSat.obs[idx_B3A];

					// 本次实验近分析B1I(北斗2代和北斗三代都有的），不考虑北三的B1(1575.42 L1(GPS))..
					//if(prn<=16)
					L1I = rnxSat.obs[idx_L1I] * t_CST::lambda(t_frequency::C2, 0);
					//if(prn>16)
					//L1I = rnxSat.obs[idx_L1I] * t_CST::lambda(t_frequency::E1, 0); ///???

					L1C = rnxSat.obs[idx_L1C] * t_CST::lambda(t_frequency::E1, 0);
					L1A = rnxSat.obs[idx_L1A] * t_CST::lambda(t_frequency::E1, 0);

					L2a = rnxSat.obs[idx_L2a] * t_CST::lambda(t_frequency::E5, 0);
					L2b = rnxSat.obs[idx_L2b] * t_CST::lambda(t_frequency::C7, 0);
					L2 = rnxSat.obs[idx_L2] * t_CST::lambda(t_frequency::E8, 0);

					L3I = rnxSat.obs[idx_L3I] * t_CST::lambda(t_frequency::C6, 0);
					L3Q = rnxSat.obs[idx_L3Q] * t_CST::lambda(t_frequency::C6, 0);
					L3A = rnxSat.obs[idx_L3A] * t_CST::lambda(t_frequency::C6, 0);

					//<判断观测值是否为零，周跳等情况
					//MP1 = value1 - (c*e+d*f)*value3 + 2*d*f*value4;
					//MP2 = value2 -  2*c*e*value3 + (c*e+d*f)*value4;

					//MP_B1I = B1I - (c*e+d*f)*L1I + 2*d*f*L2b;  //< 12 频点组合
					MP_B1I = B1I - (c4 * e4 + d4 * f4) * L1I + 2 * d4 * f4 * L3I;//<1 3频点线性组合

					//if ((abs(MP_B1I) )> 200){MP_B1I = 0; }     //txy
					//MP_B1C = B1C - (c2 * e2 + d2 * f2) * L1C + 2 * d2 * f2 * L2b;
					if (B1C != 0.0) MP_B1C = B1C - (c3 * e3 + d3 * f3) * L1C + 2 * d3 * f3 * L2a; else MP_B1C = 0;
					MP_B1A = B1A - (c2 * e2 + d2 * f2) * L1A + 2 * d2 * f2 * L2b;

					
					//if (B2a != 0.0) { MP_B2a = B2a - 2 * c1 * e1 * L1I + (c1 * e1 + d1 * f1) * L2a; }
					if (B2a != 0.0) { MP_B2a = B2a - 2 * c3 * e3 * L1C + (c3 * e3 + d3 * f3) * L2a; }
					else MP_B2a = 0;//BDS3
					if (B2b != 0.0) { MP_B2b = B2b - 2 * c * e * L1I + (c * e + d * f) * L2b; }        //BDS2
					else MP_B2b = 0;

					MP_B3I = B3I - 2 * c4 * e4 * L1I + (c4 * e4 + d4 * f4) * L3I;
					//if ((abs(MP_B3I)) > 200) { MP_B3I = 0; }
					MP_B3Q = B3Q - 2 * c4 * e4 * L1I + (c4 * e4 + d4 * f4) * L3Q;
					MP_B3A = B3A - 2 * c4 * e4 * L1I + (c4 * e4 + d4 * f4) * L3A;

					if (prn <= 16)
					{
						//MP_B2b = B2b - 2 * c * e * L1I + (c * e + d * f) * L2b;         //BDS2
						if (L1I == 0 || L3I == 0 || L2b == 0) { LMP = 0; }
					    else { LMP = (c * e - c4 * e4) * L1I + (d4 * f4) * L3I - (d * f) * L2b; }
					}
					else
					{
						//MP_B1I = B1I - (c6 * e6 + d7 * f7) * L1I + 2 * d7 * f7 * L3I;//<1 3频点线性组合
						//MP_B3I = B3I - 2 * c6 * e6 * L1I + (c6 * e6 + d7 * f7) * L3I;
						//MP_B2a = B2a - 2 * c1 * e1 * L1I + (c1 * e1 + d1 * f1) * L2a;   //BDS3
						//if ((abs(MP_B2a)) > 200) { MP_B2a = 0; }
						if (L1I == 0 || L3I == 0 || L2a == 0) { LMP = 0; }						
					    /*else { LMP = (c * e - c4 * e4) * L1I + (d4 * f4) * L3I - (d1 * f1) * L2a; }*/
						else { LMP = (c1 * e1 - c4 * e4) * L1I + (d4 * f4) * L3I -(d1 * f1) * L2a; }   //bds3 test
					}

					//<判断跳变，取平均		
					if (abs(MP_B1I - pre_B1I) > 1.0) { B1I_arc++; B1I_data[B1I_arc].wsize = 1; }
					double count(B1I_data[B1I_arc].wsize);
					B1I_data[B1I_arc].mean = (MP_B1I + (count - 1) * B1I_data[B1I_arc].mean) / count;

					if (abs(MP_B1C - pre_B1C) > 1.0) { B1C_arc++; B1C_data[B1C_arc].wsize = 1; }
					count = B1C_data[B1C_arc].wsize;
					B1C_data[B1C_arc].mean = (MP_B1C + (count - 1) * B1C_data[B1C_arc].mean) / count;

					if (abs(MP_B1A - pre_B1A) > 1.0) { B1A_arc++; B1A_data[B1A_arc].wsize = 1; }
					count = B1A_data[B1A_arc].wsize;
					B1A_data[B1A_arc].mean = (MP_B1A + (count - 1) * B1A_data[B1A_arc].mean) / count;

					if (abs(MP_B2a - pre_B2a) > 1.0) { B2a_arc++; B2a_data[B2a_arc].wsize = 1; }
					count = B2a_data[B2a_arc].wsize;
					B2a_data[B2a_arc].mean = (MP_B2a + (count - 1) * B2a_data[B2a_arc].mean) / count;

					if (abs(MP_B2b - pre_B2b) > 1.0) { B2b_arc++; B2b_data[B2b_arc].wsize = 1; }
					count = B2b_data[B2b_arc].wsize;
					B2b_data[B2b_arc].mean = (MP_B2b + (count - 1) * B2b_data[B2b_arc].mean) / count;

					if (abs(MP_B3I - pre_B3I) > 1.0) { B3I_arc++; B3I_data[B3I_arc].wsize = 1; }
					count = B3I_data[B3I_arc].wsize;
					B3I_data[B3I_arc].mean = (MP_B3I + (count - 1) * B3I_data[B3I_arc].mean) / count;

					if (abs(MP_B3A - pre_B3A) > 1.0) { B3A_arc++; B3A_data[B3A_arc].wsize = 1; }
					count = B3A_data[B3A_arc].wsize;
					B3A_data[B3A_arc].mean = (MP_B3A + (count - 1) * B3A_data[B3A_arc].mean) / count;

					if (abs(MP_B3Q - pre_B3Q) > 1.0) { B3Q_arc++; B3Q_data[B3Q_arc].wsize = 1; }
					count = B3Q_data[B3Q_arc].wsize;
					B3Q_data[B3Q_arc].mean = (MP_B3Q + (count - 1) * B3Q_data[B3Q_arc].mean) / count;

					if (abs(LMP - pre_LMP) > 0.1) { LMP_arc++; LMP_data[LMP_arc].wsize = 1; }
					count = LMP_data[LMP_arc].wsize;
					LMP_data[LMP_arc].mean = (LMP + (count - 1) * LMP_data[LMP_arc].mean) / count;

					pre_B1I = MP_B1I;//<更新
					pre_B1C = MP_B1C;
					pre_B1A = MP_B1A;

					pre_B2a = MP_B2a;
					pre_B2b = MP_B2b;

					pre_B3I = MP_B3I;
					pre_B3Q = MP_B3Q;
					pre_B3A = MP_B3A;
					pre_LMP = LMP;

					B1I_data[B1I_arc].wsize++;
					B1C_data[B1C_arc].wsize++;
					B1A_data[B1A_arc].wsize++;

					B2a_data[B2a_arc].wsize++;
					B2b_data[B2b_arc].wsize++;


					B3I_data[B3I_arc].wsize++;
					B3Q_data[B3Q_arc].wsize++;
					B3A_data[B3A_arc].wsize++;
					LMP_data[LMP_arc].wsize++;
					//outfile1 << prn << endl;
				}

			}

		}//end of while

		B1I_arc = 0; B1C_arc = 0; B1A_arc = 0; B2a_arc = 0; B2b_arc = 0; B3I_arc = 0; B3Q_arc = 0; B3A_arc = 0;
		pre_B1I = 0; pre_B1C = 0; pre_B1A = 0; pre_B2a = 0; pre_B2b = 0; pre_B2 = 0; pre_B3I = 0; pre_B3Q = 0; pre_B3A = 0;
		LMP_arc = 0; pre_LMP = 0;
		double dL4(.0);
		//<post process every epoch

		while (rnxEpo = back_grin.nextEpoch()) {

			curtime = rnxEpo->tt;//<get current obs time

			for (int i = 0; i < rnxEpo->rnxSat.size(); i++) {//<遍历各颗卫星，取得相应观测量			
				const t_rnxSat& rnxSat = rnxEpo->rnxSat[i];
				/*prn = rnxSat.satNum;*/
				if (rnxSat.satSys == satsys && rnxSat.satNum == prn) {

					B1I = rnxSat.obs[idx_B1I];
					if (B1I == 0) continue;

					B1C = rnxSat.obs[idx_B1C];
					B1A = rnxSat.obs[idx_B1A];

					B2a = rnxSat.obs[idx_B2a];
					B2b = rnxSat.obs[idx_B2b];
					B2 = rnxSat.obs[idx_B2];

					B3I = rnxSat.obs[idx_B3I];
					if (B3I == 0) continue;

					B3Q = rnxSat.obs[idx_B3Q];
					B3A = rnxSat.obs[idx_B3A];

					/*同上，不分析北三的B1*/
					L1I = rnxSat.obs[idx_L1I] * t_CST::lambda(t_frequency::C2, 0);
					//if (prn>16)  
					//L1I = rnxSat.obs[idx_L1I] * t_CST::lambda(t_frequency::E1, 0);
					L1C = rnxSat.obs[idx_L1C] * t_CST::lambda(t_frequency::E1, 0);
					L1A = rnxSat.obs[idx_L1A] * t_CST::lambda(t_frequency::E1, 0);

					L2a = rnxSat.obs[idx_L2a] * t_CST::lambda(t_frequency::E5, 0);
					L2b = rnxSat.obs[idx_L2b] * t_CST::lambda(t_frequency::C7, 0);
					L2 = rnxSat.obs[idx_L2] * t_CST::lambda(t_frequency::E8, 0);

					L3I = rnxSat.obs[idx_L3I] * t_CST::lambda(t_frequency::C6, 0);
					L3Q = rnxSat.obs[idx_L3Q] * t_CST::lambda(t_frequency::C6, 0);
					L3A = rnxSat.obs[idx_L3A] * t_CST::lambda(t_frequency::C6, 0);

					S1I = rnxSat.obs[idx_S1I];
					S1C = rnxSat.obs[idx_S1C];
					S1A = rnxSat.obs[idx_S1A];

					S2a = rnxSat.obs[idx_S2a];
					S2b = rnxSat.obs[idx_S2b];
					//S2 = rnxSat.obs[idx_B2];

					S3I = rnxSat.obs[idx_S3I];
					S3Q = rnxSat.obs[idx_S3Q];
					S3A = rnxSat.obs[idx_S3A];

					//<判断观测值是否为零，周跳等情况
					//MP1 = value1 - (c*e+d*f)*value3 + 2*d*f*value4;
					//MP2 = value2 -  2*c*e*value3 + (c*e+d*f)*value4;

					//MP_B1I = B1I - (c*e+d*f)*L1I + 2*d*f*L2b;
					MP_B1I = B1I - (c4 * e4 + d4 * f4) * L1I + 2 * d4 * f4 * L3I;//<1 3频点线性组合
					//if ((abs(MP_B1I) )> 200){MP_B1I = 0; }     //txy
					//MP_B1C = B1C - (c2 * e2 + d2 * f2) * L1C + 2 * d2 * f2 * L2b;
					if (B1C != 0) MP_B1C = B1C - (c3 * e3 + d3 * f3) * L1C + 2 * d3 * f3 * L2a; else MP_B1C = 0;
					MP_B1A = B1A - (c2 * e2 + d2 * f2) * L1A + 2 * d2 * f2 * L2b;

					//if (B2a != 0.0) MP_B2a = B2a - 2 * c1 * e1 * L1I + (c1 * e1 + d1 * f1) * L2a; //BDS3 
					if (B2a != 0.0)  MP_B2a = B2a - 2 * c3 * e3 * L1C + (c3 * e3 + d3 * f3) * L2a; 
					else MP_B2a = 0;
					if (B2b != 0.0)  MP_B2b = B2b - 2 * c * e * L1I + (c * e + d * f) * L2b;         //BDS2
					else MP_B2b=0;

					MP_B3I = B3I - 2 * c4 * e4 * L1I + (c4 * e4 + d4 * f4) * L3I;
					//if ((abs(MP_B3I)) > 200) { MP_B3I = 0; }
					MP_B3Q = B3Q - 2 * c4 * e4 * L1I + (c4 * e4 + d4 * f4) * L3Q;
					MP_B3A = B3A - 2 * c4 * e4 * L1I + (c4 * e4 + d4 * f4) * L3A;

					if (prn <= 16)
					{
						//MP_B2a = 0;   //BDS3
						//MP_B2b = B2b - 2 * c * e * L1I + (c * e + d * f) * L2b;         //BDS2
						if (L1I == 0 || L3I == 0 || L2b == 0) { LMP = 0; }
						else { LMP = (c * e - c4 * e4) * L1I + (d4 * f4) * L3I - (d * f) * L2b; }
					}
					else 
					{
						
						//MP_B1I = B1I - (c6 * e6 + d7 * f7) * L1I + 2 * d7 * f7 * L3I;//<1 3频点线性组合
						//MP_B3I = B3I - 2 * c6 * e6 * L1I + (c6 * e6 + d7 * f7) * L3I;
						//MP_B2a = B2a - 2 * c1 * e1 * L1I + (c1 * e1 + d1 * f1) * L2a;   //BDS3
						//MP_B2b = 0;//BDS2
						if (L1I == 0 || L3I == 0 || L2a == 0) { LMP = 0; }
						/*else { LMP = (c * e - c4 * e4) * L1I + (d4 * f4) * L3I - (d1 * f1) * L2a; } //L1I的频率变1575.42*/
						else { LMP = (c1 * e1 - c4 * e4) * L1I + (d4 * f4) * L3I -(d1 * f1) * L2a; }   //bds3 test
					}

					//<判断跳变，取平均		
					if (abs(MP_B1I - pre_B1I) > 1.0) B1I_arc++;
					if (abs(MP_B1C - pre_B1C) > 1.0) B1C_arc++;
					if (abs(MP_B1A - pre_B1A) > 1.0) B1A_arc++;
					if (abs(MP_B2a - pre_B2a) > 1.0) B2a_arc++;
					if (abs(MP_B2b - pre_B2b) > 1.0) B2b_arc++;
					if (abs(MP_B3I - pre_B3I) > 1.0) B3I_arc++;
					if (abs(MP_B3Q - pre_B3Q) > 1.0) B3Q_arc++;
					if (abs(MP_B3A - pre_B3A) > 1.0) B3A_arc++;
					if (abs(LMP - pre_LMP) > 0.1) LMP_arc++;   //txy

					//if (abs(dL1IL2a - pre_dL1IL2a) > 0.05) { L1I_arc++; L2a_arc++; }
					////if (abs(dL1IL2b -pre_dL1IL2b)>0.15) {L1I_arc++; L2b_arc++;} 
					//dL4 = dL1IL2a - pre_dL1IL2a;
					//pre_dL1IL2a = dL1IL2a;
					//pre_dL1IL2b = dL1IL2b;



					pre_B1I = MP_B1I;//<更新
					pre_B1C = MP_B1C;
					pre_B1A = MP_B1A;

					pre_B2a = MP_B2a;
					pre_B2b = MP_B2b;

					pre_B3I = MP_B3I;
					pre_B3Q = MP_B3Q;
					pre_B3A = MP_B3A;
					pre_LMP = LMP;  //txy GFIF组合



					//<2------------get satellite position----------------------------
					double xc[4], vv[3];
					ColumnVector satpos(3);
					bncTime bt(curtime - (B3I) / t_CST::c), tt(curtime - (B3I) / t_CST::c);//signal broadcast time

					bdseph = rnxeph.getBDSEph(bt, rnxSat.satNum);
					bdseph.position(bt.gpsw(), bt.gpssec(), xc, vv);

					for (int i = 0; i < 1; i++) {
						bt = tt - xc[3]; //< clk+rel 卫星钟差改正			
						bdseph.position(bt.gpsw(), bt.gpssec(), xc, vv);//<读取星历时已进行时间转换
					}

					satpos[0] = xc[0]; satpos[1] = xc[1]; satpos[2] = xc[2];
					double elevation(0.0), azimuth(0.0);
					cmpEle(TruePos, satpos, elevation, azimuth);



					//outfile<<curtime.datestr(' ')<<"  "<<curtime.timestr(0,' ')<<"    "<<prn
					//	<<"   "<<elevation/t_CST::M_PI*180
					//	<<"   "<<B1I_arc<<"   "<<MP_B1I-B1I_data[B1I_arc].mean
					//	<<"   "<<B1C_arc<<"   "<<MP_B1C-B1C_data[B1C_arc].mean
					//	<<"   "<<B1A_arc<<"   "<<MP_B1A-B1A_data[B1A_arc].mean
					//	<<"   "<<B2a_arc<<"   "/*<<L2a_arc<<"    "<<dL4<<"   "*/<<MP_B2a-B2a_data[B2a_arc].mean
					//	<<"   "<<B2b_arc<<"   "<<MP_B2b-B2b_data[B2b_arc].mean
					//	<<"   "<<B3I_arc<<"   "<<MP_B3I-B3I_data[B3I_arc].mean
					//	<<"   "<<B3Q_arc<<"   "<<MP_B3Q-B3Q_data[B3Q_arc].mean
					//	<<"   "<<B3A_arc<<"   "<<MP_B3A-B3A_data[B3A_arc].mean
					//	<<"  "<<S1I<<"  "<<S1C<<"  "<<S1A<<"  "<<S2a
					//	<<"  "<<S2b<<"  "<<S3I<<"  "<<S3Q<<"  "<<S3A<<"   "<<endl;

				//	outfile << curtime.datestr(' ') << "  " << curtime.timestr(0, ' ') << "    " << prn
				//		<< setiosflags(ios::left) << setw(8) << elevation / t_CST::M_PI * 180;
				//	if (MP_B1I) {
				//		outfile << setiosflags(ios::left) << setw(8) << MP_B1I - B1I_data[B1I_arc].mean
				//			<< setiosflags(ios::left) << setw(8) << S1I;
				//	}
				//	if (MP_B1C)
				//	{
				//		outfile << setiosflags(ios::left) << setw(8) << MP_B1C - B1C_data[B1C_arc].mean<< setiosflags(ios::left) << setw(8) << S1C;
				//}
				//	if (MP_B1A)
				//	{
				//		outfile << setiosflags(ios::left) << setw(8) << MP_B1A - B1A_data[B1A_arc].mean<< setiosflags(ios::left) << setw(8) << S1A;
				//	}
				//	if (MP_B2a) {
				//		outfile << setiosflags(ios::left) << setw(8) << MP_B2a - B2a_data[B2a_arc].mean << setiosflags(ios::left) << setw(8) << S2a;
				//	}
				//	if (MP_B2b) {
				//		outfile << setiosflags(ios::left) << setw(8) << MP_B2b - B2b_data[B2b_arc].mean << setiosflags(ios::left) << setw(8) << S2b;
				//	}
				//	if (MP_B3I) {
				//		outfile << setiosflags(ios::left) << setw(8) << MP_B3I - B3I_data[B3I_arc].mean << setiosflags(ios::left) << setw(8) << S3I;
				//	}

				//	if (MP_B3Q) {
				//		outfile << setiosflags(ios::left) << setw(8) << MP_B3Q - B3Q_data[B3Q_arc].mean << setiosflags(ios::left) << setw(14) << S3Q;
				//	}
				//	if (MP_B3A) {
				//		outfile << setiosflags(ios::left) << setw(8) << MP_B3A - B3A_data[B3A_arc].mean << setiosflags(ios::left) << setw(14) << S3A;
				//	}
				//	outfile << endl;
				/****           output C01 form      *****/
					//string c("0"), prns;
					//if (prn < 10)
					//{						
					//	prns = to_string(prn);
					//	prns = c+prns;
					//}
					//else 
					//{ prns = to_string(prn); }
			   /****           output C01 form        ******/
					 
					     //----------------        rnx304 bds 有C1 C6     ---------------------- //
						outfile << curtime.datestr(' ') << "  " << curtime.timestr(0, ' ') << "    " << setiosflags(ios::right) << setw(2) <<prn//<< "C" << prns << "    "
							<< setiosflags(ios::right) << setprecision(3) << setw(10) << azimuth / t_CST::M_PI * 180
							<< setiosflags(ios::right) << setprecision(3) << setw(10) << elevation / t_CST::M_PI * 180
							<< setiosflags(ios::right) << setw(8) << setprecision(2) << S1I
							<< setiosflags(ios::right) << setw(8) << setprecision(3) << MP_B1I - B1I_data[B1I_arc].mean
							<< setiosflags(ios::right) << setw(8) << setprecision(2) << S1C
							<< setiosflags(ios::right) << setw(8) << setprecision(3) << MP_B1C - B1C_data[B1C_arc].mean
							//<< setiosflags(ios::right) << setw(8) <<MP_B1C-B1C_data[B1C_arc].mean
							//<< setiosflags(ios::right) << setw(8) <<MP_B1A-B1A_data[B1A_arc].mean 
							<< setiosflags(ios::right) << setw(8) << setprecision(2) << S2a
							<< setiosflags(ios::right) << setprecision(3) << setw(8) << MP_B2a - B2a_data[B2a_arc].mean
							<< setiosflags(ios::right) << setw(8) << setprecision(2) << S2b
							<< setiosflags(ios::right) << setprecision(3) << setw(8) << MP_B2b - B2b_data[B2b_arc].mean
							<< setiosflags(ios::right) << setprecision(2) << setw(8) << S3I
							<< setiosflags(ios::right) << setprecision(3) << setw(8) << MP_B3I - B3I_data[B3I_arc].mean
							<< setiosflags(ios::right) << setprecision(3) << setw(8) << LMP - LMP_data[LMP_arc].mean << endl;
						 //---------------      rnx302 bds 只有C1 C7 不含C6    --------------------//
						//outfile << curtime.datestr(' ') << "  " << curtime.timestr(0, ' ') << "    " << setiosflags(ios::right) << setw(1) << "C" << prns << "    "
						//	<< setiosflags(ios::right) << setprecision(3) << setw(10) << azimuth / t_CST::M_PI * 180
						//	<< setiosflags(ios::right) << setw(8) << elevation / t_CST::M_PI * 180
						//	<< setiosflags(ios::right) << setw(8) << setprecision(2) << S1I
						//	<< setiosflags(ios::right) << setw(8) << setprecision(3) << MP_B1I - B1I_data[B1I_arc].mean
						//	//<< setiosflags(ios::right) << setw(8) <<MP_B1C-B1C_data[B1C_arc].mean
						//	//<< setiosflags(ios::right) << setw(8) <<MP_B1A-B1A_data[B1A_arc].mean 
						//	<< setiosflags(ios::right) << setprecision(3) << setw(8) << ""
						//	<< setiosflags(ios::right) << setprecision(3) << setw(8) << ""
						//	<< setiosflags(ios::right) << setprecision(2) << setw(8) << S2b
						//	<< setiosflags(ios::right) << setprecision(3) << setw(8) << MP_B2b - B2b_data[B2b_arc].mean
						//	<< setiosflags(ios::right) << setprecision(3) << setw(8) << LMP - LMP_data[LMP_arc].mean << endl;
						
						//if (prn <= 16)  //BDS3  txy
					//	outfile << setiosflags(ios::right) << setprecision(2) << setw(8) << S2b << endl;
					//if (prn > 16)
					//	outfile << setiosflags(ios::right) << setprecision(2) << setw(8) << S2a << endl;
					

					//<< setiosflags(ios::right) << setw(8) <<MP_B3Q-B3Q_data[B3Q_arc].mean
					//<< setiosflags(ios::right) << setw(8) <<MP_B3A-B3A_data[B3A_arc].mean					
					//<< setiosflags(ios::right) << setw(14) << S1C
					//<< setiosflags(ios::right) << setw(14) << S1A
					//<< setiosflags(ios::right) << setw(14) << S2a
					//<< setiosflags(ios::right) << setw(14) << S2b					
					//<< setiosflags(ios::right) << setw(14) << S3Q
					//<< setiosflags(ios::right) << setw(14) << S3A << endl;
					//cout << curtime.timestr(0, ' ');
				}
			}

		
		}//end of while
		cout << prn;
	}
	return 0;
}