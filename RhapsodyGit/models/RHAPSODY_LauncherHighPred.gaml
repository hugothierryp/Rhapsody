/**
 *  RHAPSODY_Launcher
 *  Author: Hugo THIERRY
 *  Description: 
 */

model RHAPSODY_Launcher

import 'ATLAS_RHAPSODYVolunteerHighPred.gaml'
import 'RHAPSODYVolunteerHighPred.gaml'

global torus:true { 
	
	// All interface parameters are declared in the global in GAMA
	file InitialisationFile <- csv_file('../includes/VCGRHAPSODYLaunchHighPred.csv');	 //Default set of parameters to initialize the simulation
	file landscapeGIS <- shape_file('../includes/Parcellaire2012Eclater.shp');	// GIS File
	file climateCsv;						// Wind and Temperature File
	file rotationsCSV;		// Rotations File
	file landcoverCSV;					// Covers File
	file landuseCSV;					// Covers File	
	file cropStagesCSV;				//Crop Stages File											
	file rainCSV;							// Rain File
	
	int rep<-1;  /// Parameter used for the batch outputs
	
	file rasterFile <- grid_file('../includes/vcgraster30m.txt');				// GIS File input
	file immigrationData;			// Aphid Immigration data
	file aphidCropRates;			// Aphid reproduction rates depending on crop stages
	file aphidMortalityRates;
	
	file mask1<-csv_file('../includes/Mask1.csv');
	file mask2<-csv_file('../includes/Mask2.csv');
	file mask3<-csv_file('../includes/Mask3.csv');
	file mask4<-csv_file('../includes/Mask4.csv');
	file mask5<-csv_file('../includes/Mask5.csv');
	file mask6<-csv_file('../includes/Mask6.csv');
	file mask7<-csv_file('../includes/Mask7.csv');
	file mask8<-csv_file('../includes/Mask8.csv');		
	
	bool hasHourlyWeatherData;
	bool hasInputRotations; 				// Do we have already a rotation and an index of succession for each agricultural patch in the gis file ?
	
	int timestepHours; 						// Hour many hours does one timestep represent in the simulation	
	int initialDay;
	int initialMonth;
	int julianDay;								// Day of the year
	int daysElapsed;							// Number of days since the simulation started
	int hour;									// Hour of the day	
	int year;									// Year of the simulation
	
	
	Landscape myLandscape; //Create the landscape agent
	
	geometry shape <- envelope(landscapeGIS);	// Limit the model to the GIS boundaries
	
	
	init{
		matrix fichiersData <- matrix(InitialisationFile);

		climateCsv <- csv_file('../includes/'+string(fichiersData[0,0]));
		rotationsCSV <- csv_file('../includes/'+string(fichiersData[0,1]));
		landcoverCSV <- csv_file('../includes/'+string(fichiersData[0,2]));
		landuseCSV <- csv_file('../includes/'+string(fichiersData[0,3]));
		cropStagesCSV <- csv_file('../includes/'+string(fichiersData[0,4]));
		rainCSV <- csv_file('../includes/'+string(fichiersData[0,5]));
		initialDay<-int(fichiersData[0,6]);
		initialMonth<-int(fichiersData[0,7]);
		hasHourlyWeatherData<-bool(fichiersData[0,8]);
		if hasHourlyWeatherData {
			timestepHours<-int(fichiersData[0,9]);
		}
		else{
			timestepHours<-8;
		}
		hasInputRotations<-bool(fichiersData[0,10]);
		immigrationData<- csv_file('../includes/'+string(fichiersData[0,11]));	
		aphidCropRates<- csv_file('../includes/'+string(fichiersData[0,12]));
		
		list<int> daysPerMonth <- [0,31,59,90,120,151,181,212,243,273,304,334];
		julianDay<-daysPerMonth[initialMonth-1]+initialDay;

		hour <-24-timestepHours; 		//Start just before the first day
		daysElapsed <--1;
		year <-0;						 // Start at year = 0;
		
		create Landscape number:1 returns: createdLandscape; 	// Creates the landscape 
		myLandscape<- any(createdLandscape); 					//createdLandscape is a list of one entity
	}
			
	reflex update_time{
	//Purpose: Update the year, day and hour of the simulation
	// InOut :self.year, self.day, self.hour, self.numberOfDays			
		hour <- hour +timestepHours;			
		if hour = 24 {hour <- 0; julianDay <- julianDay + 1; daysElapsed <- daysElapsed +1;}
		if julianDay = 366 {year <- year+1;julianDay <- 1;}	
	}

			
}

experiment RHAPSODY_Landscape type: gui{

	parameter 'Select the initiailization file defining the simulation' var: InitialisationFile category: 'File input';
			
	output {
		display patch_display refresh: every(cycle) {										// Shows a 2D image of the landscape
			species Patch aspect: base ;
			grid raster_landscape lines: #black transparency:0.3 ;								
		}
		display myDisplay type:opengl {													// Show a 3D representation of the landscape
			 species Patch aspect: base;
		}
		display "Daily Crop Surfaces" refresh: every(cycle){
			chart "Daily Crop Surfaces" type: series background: #white {
				data "Wheat" value: Landscape(0).totalWheat color: rgb("yellow") line_visible: false;
				data "Corn" value: Landscape(0).totalCorn color: rgb("orange") line_visible: false;
				data "Volunteer" value: Landscape(0).totalWheatVolunteer color: rgb("black") line_visible: true;
				data "Sorghum" value: Landscape(0).totalSorghum color: rgb("red")line_visible: false  ;							
			}
		}		
	   display Climat refresh: every(cycle) {
			chart "Mean Temperature" type: series background: rgb("white") size: {1,0.4} position: {0, 0.05}  {			// Mean temperature of the landscape 
				data "Temperature" value: Landscape(0).climate.landscapeTemperature   color: rgb("blue") ;
			}
			chart "Wind Speed" type: series background: rgb("white") size: {1,0.4} position: {0, 1000}  {				// Mean windspeed of the landscape
				data "WindSpeed" value: Landscape(0).climate.windSpeed color: rgb("blue") ;
			}
		}	
			
		monitor julianDay value: julianDay refresh: every(1 #cycle) ;
		monitor Days_elapsed value:daysElapsed refresh: every(1 #cycle);
		monitor Hour value: hour refresh: every(1 #cycle) ;
		monitor Year value: year refresh: every(1 #cycle
		
		);

	}	
}

experiment RHAPSODY_Batch type: batch repeat:1 keep_seed: false  until: (year > 11){
	parameter 'Repetition' var: rep among: [ 1,2,3,4,5,6,7,8,9,10];
}

