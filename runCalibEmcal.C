#!/usr/bin/gawk {system("root -b -q " FILENAME); exit;}
// -*- mode: c++; -*-
//#include <TROOT.h>
//#include <TSystem.h>
#define RUN_PERIOD(x, y) ((x) << 5) + ((y) - 'a')


void runCalibEmcal(const char *run_mode = "full", const int lhc_run_period = RUN_PERIOD(18, 'm'))
/*/////////////////////////////////////////////////////////////////////////////////
Run Mode:
- "test": runs the code on a single node, and with a small set of events,
afterwards it will transfer the ROOT format output file back

- "full": runs the code entirely, with the run numbers you specified, and
afterwards it will attempt to transfer back all the ROOT format output
files and merge them (though you can do both manually). Note that "full"
will fail, if your output directory is already occupied, and you try to
submit on top of the existing directory.

- "terminate": only does the transfer of the output, and merge them.

Run period(number, letter)
*//////////////////////////////////////////////////////////////////////////////////
{
	gROOT->ProcessLine(".include $ROOTSYS/include");
	gROOT->ProcessLine(".include $ALICE_ROOT/include");
	gROOT->ProcessLine(".include $ALICE_PHYSICS/include");

	// Load base root libraries
	gSystem->Load("libTree");
	gSystem->Load("libGeom");
	gSystem->Load("libVMC");
	gSystem->Load("libSTEERBase");
	gSystem->Load("libPhysics");
	gSystem->Load("libESD");
	gSystem->Load("libAOD");

	// Load analysis framework libraries
	gSystem->Load("libANALYSIS");
	gSystem->Load("libANALYSISalice");
	gSystem->Load("libEMCALUtils");
	gSystem->Load("libPWGPPEMCAL");

	gROOT->ProcessLine(".L AliAnalysisTaskCalibEmcal.cxx+g");

	AliAnalysisManager *mgr = new AliAnalysisManager();

	gROOT->Macro("$ALICE_ROOT/ANALYSIS/macros/train/"
				 "AddESDHandler.C");

	AliAnalysisAlien *plugin =
		new AliAnalysisAlien("pluginCalibEmcal");

	plugin->AddIncludePath("-I. -I$ROOTSYS/include -I$ALICE_ROOT "
			       "-I$ALICE_ROOT/include -I$ALICE_PHYSICS/include");

	plugin->SetGridWorkingDir("workdir_calib");
	plugin->SetGridOutputDir("outputdir_18m_14runsPart14");
	plugin->SetAliPhysicsVersion("vAN-20190109-1");
	plugin->SetAdditionalLibs(
		"AliAnalysisTaskCalibEmcal.h "
		"AliAnalysisTaskCalibEmcal.cxx");
	plugin->SetAnalysisSource("AliAnalysisTaskCalibEmcal.cxx");
	plugin->SetRunPrefix("000");

	const int run_number_lhc10h[] = {

		137135, 137161, 137162, 137230, 137231, 137232, 137235,
		137236, 137243, 137366, 137430, 137431, 137432, 137434,
		137439, 137440, 137441, 137443, 137530, 137531, 137539,
		137541, 137544, 137546, 137549, 137595, 137608, 137638,
		137639, 137685, 137686, 137691, 137692, 137693, 137704,
		137718, 137722, 137724, 137751, 137752, 137844, 137848,
		138190, 138192, 138197, 138201, 138225, 138275, 138364,
		138396, 138438, 138439, 138442, 138469, 138534, 138578,
		138582, 138583, 138621, 138624, 138638, 138652, 138653,
		138662, 138666, 138730, 138732, 138837, 138870, 138871,
		138872, 139028, 139029, 139036, 139037, 139038, 139105,
		139107, 139173, 139309, 139310, 139314, 139328, 139329,
		139360, 139437, 139438, 139465, 139503, 139505, 139507,
		139510,

		-1
	};

	const int run_number_lhc11h[] = {

		170593, 170572, 170388, 170312, 170311, 170309, 170308,
#if 1
		170270, 170269, 170268, 170230, 170228, 170207, 170204,
		170203, 170193, 170155, 170083, 169965, 169923, 169859,
		169858, 169855, 169846, 169838, 169837, 169835, 169591,
		169588, 169557, 169554, 169550, 169512, 169504, 169498,
		169475, 169419, 169417, 169415, 169411, 169238, 169160,
		169156, 169144, 169138, 169035, 168826, 168512, 168511,
		168467, 168464, 168458, 168361, 167988, 167987, 167920,
		167915,
#endif

		-1
	};

	const int run_number_lhc12b[] = {

		177580, 177592, 177597, 177601, 177624, 177681, 177798,
		177799, 177804, 177805, 177810, 177858, 177860, 177861,
		177864, 177869, 177932, 177938, 177942, 178018, 178024,
		178025, 178026, 178028, 178029, 178030, 178031, 178052,
		178053, 178163, 178167, 178220,

		-1
	};

	const int run_number_lhc15f[] = {

		224895, 224896, 224897, 224898, 224930, 224988, 224997,
		225000, 225011, 225026, 225031, 225035, 225050, 225051,
		225052, 225093, 225105, 225106, 225576, 225578, 225579,
		225580, 225582, 225586, 225587, 225702, 225709, 225716,
		225763, 225768, 225799, 225917, 226114, 226115, 226166,
		226440, 226444, 226445, 226466, 226468, 226530,
		226532,

		-1
	};

	const int run_number_lhc15i[] = {

		235716, 235812, 235840, 235887, 235894, 235899, 236149,
		236282, 236332, 236350, 236355, 236358, 236141, 236156,
		236162, 236225, 236232, 236241, 236245, 236391, 236394,
		236442, 236445, 236447, 236454, 236466, 236555, 236559,
		236819, 236820, 236823, 236849,

		-1
	};
	const int run_number_lhc15j[] = {

		236967, 236970, 237003, 237030, 237050, 237106, 237109,
		237119, 237178, 237257, 237287, 237511, 237513, 237516,
		237646, 237673, 237686, 237700, 237701, 237703, 237704,
		237709, 237712, 237766, 237773, 237778, 237781, 237788,
		237792, 237794, 237805, 237843, 237846,

		-1
	};

	const int run_number_lhc15o[] = {
	  
	  //Globally good runs without those with special conditions, MB, or special bad channel map
	  /*245146, 245151, 245152, 245231, 245232, 245259, 245343, 
	  245345, 245346, 245347, 245349, 245353, 245396, 245397, 
	  245401, 245407, 245409, 245411, 245439, 245441, 245446, 
	  245454, 245496, 245497, 245501, 245504, 245505, 245507, 
	  245535, 245540, 245542, 245543, 245544, 245545, 245554, 
	  245700, 245702, 245705, 245738, 245829, 245831, 245833, 
	  245949, 245952, 245954, 245963, 246001, 246003, 246037, 
	  246042, 246052, 246053, 246087, 246089, 246113, 246115, 
	  246217, 246222, 246225, 246271, 246272, 246390, 246391, 
	  246392, 246424, 246434, 246487, 246488, 246493, 246495, 
	  246750, 246751, 246757, 246758, 246759, 246760, 246765, 
	  246766, 246804, 246805, 246807, 246808, 246809, 246810, 
	  246844, 246845, 246846, 246928, 246945,//*/
	  
	  //List 1
	  /*244918 , 244975 , 244980 , 244982 , 244983 , 245064 , 
	  245066 , 245068 , 245145 , 245146 , 245151 , 245152 , 
	  245231 , 245232 , 245259 , 245343 , 245345 , 245346 , 
	  245347 , 245349 , 245353 , 245396 , 245397 , 245401 , 
	  245407 , 245409 , 245411 , 245439 , 245441 , 245446 ,//new list*/

	  /*245683, 245700, 245702, 245705, 245829, 245831, 245833,
	  245949, 245952, 245954, 246001, 246003, 246037, 246042,
	  246052, 246053, 246087, 246089, 246113, 246115,//*/
	  
	  
	  //List 2
	  245496, 245497, 245501, 245504, 245505, 245507, 245535, 245540, 245542, 245543, 245544,245545 , 245554, 245683, 245700, 245702, 245705, 245738, 245829, 245831, 245833, 245949, 245952, 245954, 245963, 246001, 246003, 246037, 246042, 246052, 246053, 246087, 245996, 245729, 245731, 245785,//*/ new list
	  //246217, 246222, 246225, 246271, 246272,
	  

	  //List 3
	  /*245554 , 245683 , 245700 , 245702 , 245705 , 245738 , 
	  245829 , 245831 , 245833 , 245949 , 245952 , 245954 , 
	  245963 , 246001 , 246003 , 246037 , 246042 , 246052 , 
	  246053 , 246087 , 246089 , 246113 , 246115 , 246217 , 
	  246222 , 246225 , 246271 ,//new list*/
	  
	  /*246242, 246433, 246434, 246487, 246488, 246493, 246495,
	  246540, 246543, 246567, 246568, 246575, 246583, 246648,
	  246750, 246751, 246757, 246758, 246759, 246760, 246765,
	  246766, 246804, 246805, 246807, 246808, 246809, 246810,
	  246844, 246845, 246846,//*/
	  

	  //list 5
	  //246272 , 246390 , 246391 , 246392 , 246424 , 246434 , 246487 , 246488 , 246493 , 246495 , 246750 , 246751 , 246757 , 246758 , 246759 , 246760 , 246765 , 246766 , 246804 , 246805 , 246807 , 246808 , 246809 , 246810 , 246844 , 246845 , 246846 , 246928 , 246945 ,//new list
	  //246928, 246930, 246945,
	  

		-1
	};

	const int run_number_lhc16h[] = {
#if 1
		// Restricted list
		254604, 254606, 254607, 254621, 254629, 254630, 254632,
		254640, 254644, 254646, 254648, 254649, 254651, 254652,
		254653, 254654,
#else
		254604, 254606, 254607, 254621, 254629, 254630, 254632,
		254640, 254644, 254646, 254648, 254649, 254651, 254652,
		254653, 254654, 255249, 255251, 255252, 255253, 255255,
		255256, 255275, 255276, 255350, 255351, 255352, 255418,
		255419, 255420, 255421, 255440, 255463, 255465, 255466,
		255467,
#endif

		-1
	};
	
	const int run_number_lhc16q[] = {
	  //Globally Good
	  265309, 265332, 265334, 265335, 265336, 265338,
	  265339, 265342, 265343, 265344, 265378, 265383, 
	  265384, 265387, 265388, 265419, 265420, 265421, 
	  265422, 265424, 265425, 265426, 265427, 265435, 
	  265499, 265500, 265501, 265521, 265525,
	  

	  -1
	};


	const int run_number_lhc16r[] = {
	  //Globally Good
	  266318, 266317, 266316, 266208, 266197, 266196, 
	  266187, 265744,
	  
	  -1
	};
	
	const int run_number_lhc16s[] = {
	  //Globally Good
	  267110, 267081, 267077, 267072, 267070, 266998, 
	  266997, 266994, 266993, 266944, 266886, 266885, 
	  266883, 266882, 266437,

	  -1
	};

	const int run_number_lhc16t[] = {
	  //Globally Good
	  267163, 267164, 267165, 267166,

	  -1
	  };

	
	const int run_number_lhc17p[] = {
	  
	  282343, 282342, 282341, 282340, 282314, 282313, 282312,
	  282307, 282306, 282305, 282304, 282303, 282302, 282247,
	  282230, 282229, 282227, 282224, 282206, 282189, 282147,
	  282146, 282127, 282126, 282125, 282123, 282122, 282119,
	  282118, 282099, 282098, 282078, 282051, 282031, 282030,
	  282025,
	  
		-1
	};

	const int run_number_lhc18d[] = {
	  
	  286014, 286018, 286025, 286026, 286027, 286030, 286064,
	  286124, 286127, 286129, 286130, 286159, 286198, 286201,
	  286202, 286203, 286229, 286230, 286231, 286254, 286255,
	  286256, 286257, 286258, 286261, 286263, 286282, 286308,
	  286309, 286310, 286311, 286314, 286348, 286349,
	  
	  -1
	};

	const int run_number_lhc18m[] = {
	  
	  /*290293, 290294, 290297,
	  290298, 290300, 290323, 290324, 290327, 290350, 290374,
	  290375, 290376, 290399, 290401, 290411, 290412, 290418,
	  290420, 290421, 290425, 290426, 290427, 290428, 290456,
	  290458, 290459, 290469, 290499, 290500, 290533, 290534,
	  290535, 290538, 290539, 290540, 290544, 290549, 290550,
	  290553, 290588, 290590, 290612, 290613, 290614, 290615,
	  290627, 290632, 290645, 290658, 290660, 290665, 290687,
	  290689, 290692, 290696, 290699, 290721, 290742, 290764,
	  290766, 290769, 290772, 290774, 290787, 290790, 290841,
	  290843, 290846, 290848, 290860, 290862, 290886, 290887,
	  290892, 290894, 290895, 290932, 290935, 290941, 290943,
	  290944, 290948, 290974, 290975, 290976, 290979, 290980,
	  291002, 291003, 291004, 291005, 291035, 291037, 291038,
	  291041, 291065, 291066, 291069, 291093, 291100, 291101,
	  291110, 291111, 291116, 291143, 291188, 291209, 291238,
	  291240, 291257, 291262, 291263, 291265, 291266, 291282,
	  291283, 291284, 291285, 291286, 291360, 291361, 291362,
	  291373, 291375, 291377, 291397, 291399, 291402, 291417,
	  291419, 291420, 291424, 291446, 291447, 291453, 291456,
	  291457, 291481, 291482, 291484, 291485, 291613, 291614,
	  291615, 291618, 291622, 291624, 291625, 291626, 291657
	  291661, 291665, 291690, 291692, 291694, 291698, 291706,
	  291729, 291755, 291756, 291760, 291768, 291795, 291796,
	  291803, 291805, 291807, 291942, 291943, 291944, 291945,
	  291946, 291948, 291953, 291976, 291977, 291982, 292012,
	  292040, 292060, 292061, 292062, 292067, 292075, 292080,
	  292081, 292106, 292107, 292108, 292109, 292114, 292115,
	  292140, 292160, 292161, 292162, 292163, 292164, 292166,
	  292167, 292168, 292192, 292218, 292240, 292241, 292242,
	  292265, 292270, 292273, 292274, 292298, 292397, 292398,
	  292405, 292406, 292428, 292429, 292430, 292432, 292434,//*/
	  292456, 292457, 292460, 292461, 292495, 292496, 292497,
	  292500, 292521, 292523, 292524, 292526, 292553, 292554,
	  /*292557, 292559, 292560, 292563, 292584, 292586, 292593,
	  292598, 292693, 292695, 292696, 292698, 292701, 292704,
	  292737, 292739, 292747, 292748, 292750, 292809, 292810,
	  292811, 292831, 292832, 292836,//*/
	  
	  -1
	};

	const int *run_number;

	plugin->SetGridDataDir(Form(
		"/alice/data/%d/LHC%02d%c", 2000 + (lhc_run_period >> 5),
		lhc_run_period >> 5, (lhc_run_period & 0x1f) + 'a'));

	switch (lhc_run_period) {
	case RUN_PERIOD(10, 'h'):
		plugin->SetDataPattern("ESDs/pass2/*/AliESDs.root");
		run_number = run_number_lhc10h;
		break;
	case RUN_PERIOD(11, 'h'):
		plugin->SetDataPattern("ESDs/pass1_calo/*/AliESDs.root");
		run_number = run_number_lhc11h;
		break;
	case RUN_PERIOD(12, 'b'):
		plugin->SetDataPattern("pass1/*/AliESDs.root");
		run_number = run_number_lhc12b;
		break;
	case RUN_PERIOD(15, 'f'):
		plugin->SetDataPattern("pass1/*/AliESDs.root");
		run_number = run_number_lhc15f;
		break;
	case RUN_PERIOD(15, 'i'):
		plugin->SetDataPattern("muon_calo_pass2/*/AliESDs.root");
		run_number = run_number_lhc15i;
		break;
	case RUN_PERIOD(15, 'j'):
		plugin->SetDataPattern("muon_calo_pass2/*/AliESDs.root");
		run_number = run_number_lhc15j;
		break;
	case RUN_PERIOD(15, 'o'):
		plugin->SetDataPattern("muon_calo_pass1/*/AliESDs.root");
		run_number = run_number_lhc15o;
		break;
	case RUN_PERIOD(16, 'h'):
		plugin->SetDataPattern("muon_calo_pass1/*/AliESDs.root");
		run_number = run_number_lhc16h;
		break;
	case RUN_PERIOD(16, 'q'):
	        plugin->SetDataPattern("muon_calo_pass1/*/AliESDs.root");
	        run_number = run_number_lhc16q;
	        break;

	case RUN_PERIOD(16, 'r'):
	        plugin->SetDataPattern("muon_calo_pass1/*/AliESDs.root");
	        run_number = run_number_lhc16r;
	        break;

	case RUN_PERIOD(16, 's'):
	        plugin->SetDataPattern("muon_calo_pass1/*/AliESDs.root");
	        run_number = run_number_lhc16s;
	        break;

	case RUN_PERIOD(16, 't'):
	        plugin->SetDataPattern("muon_calo_pass1/*/AliESDs.root");
	        run_number = run_number_lhc16t;
	        break;
		
	case RUN_PERIOD(17, 'p'):
	        plugin->SetDataPattern("muon_calo_pass1/*/AliESDs.root");
		run_number = run_number_lhc17p;
		break;
	case RUN_PERIOD(18, 'd'):
	        plugin->SetDataPattern("muon_calo_pass1/*/AliESDs.root");
		run_number = run_number_lhc18d;
		break;
	case RUN_PERIOD(18, 'm'):
	        plugin->SetDataPattern("muon_calo_pass1/*/AliESDs.root");
		run_number = run_number_lhc18m;
		break;

	}
	




	for (const int *r = run_number; *r != -1; r++) {
		plugin->AddRunNumber(*r);
	}


	const char *alien_close_se = gSystem->Getenv("alien_CLOSE_SE");

	if (alien_close_se != NULL) {
		const char *file = mgr->GetCommonFileName();

		plugin->SetDefaultOutputs(kFALSE);
		plugin->SetOutputFiles(Form(
			"%s@%s", file, alien_close_se));
		plugin->SetOutputArchive(Form(
			"log_archive.zip:stdout,stderr@%s "
			"root_archive.zip:%s,*.stat@%s",
			alien_close_se, file, alien_close_se));
	}

	plugin->SetRunMode(run_mode);
	mgr->SetGridHandler(plugin);
	gROOT->Macro("macros/AddAliAnalysisTaskCalibEmcal.C");
	if (mgr->InitAnalysis()) {
		mgr->StartAnalysis("grid");
	}
}
