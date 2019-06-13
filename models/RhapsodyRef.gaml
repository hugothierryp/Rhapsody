/**
 *  RHAPSODY_Launcher
 *  Author: Hugo THIERRY
 *  Description: 
 */

model RHAPSODY_Launcher

global torus:true { 
	
	// All interface parameters are declared in the global in GAMA
	file InitialisationFile <- csv_file('../includes/VCGRHAPSODYLaunchRef.csv');	 //Default set of parameters to initialize the simulation
	file landscapeGIS  <- shape_file('../includes/Parcellaire2012Eclater.shp');	// GIS File
	file climateCsv;						// Wind and Temperature File
	file rotationsCSV;		// Rotations File
	file landcoverCSV;					// Covers File
	file landuseCSV;					// Covers File	
	file cropStagesCSV;				//Crop Stages File											
	file rainCSV;							// Rain File
	
	int rep<-1;
	
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

species Landscape{
	list<Landuse> landusesList<-[]; // List of the landuses
	list<Landcover> landcoversList<-[]; // List of the landcovers
	list<Rotation> rotationsList<-[]; // List of the rotations
	list<raster_landscape> availableHabitats<-[];
	list<Patch> patchesList<-[]; // List of the patches
	Climate climate; // Climate applied to the landscape
	
	int aphidImmigrationIndex;				// Index used in the aphid immigration file		
	list<float> aphidMovementDistProbs<-list_with(number_movement_circles-1,0.0);
	list<float> distanceGroupProbabilities<-list_with(number_movement_circles-1,0.0);
	matrix dailyImmigrantData <-matrix(immigrationData);			// Data matrix
	int landscape_grid_size_X <- 84;
	int landscape_grid_size_Y <- 84;
	int number_movement_circles <- 8;
	int grid_cell_size <- 30;
	int numberMovementCells <-15;
	
	
	map<int,matrix> movementMask;
	
	float totalWheat;
	float totalCorn;
	float totalWheatVolunteer;
	float totalSorghum;
			
	
	init{
	// Purpose: create the components of a Landscape from the global CSV and SHP files
	// out: self.landusesList, self.rotationsList, self.patchesList, self.climate			
		create Climate number:1 returns:createdClimate; // Create the climate entitie
		climate<-any(createdClimate); // createdClimate is a list of one entity
		do create_landuses();
		do create_landcovers();
		do create_rotations();
		do create_patches();
		ask raster_landscape{
			do define_reference_patch();
			do define_edge_cells();			
		}
		do calculate_movement_probabilities();
		ask patchesList where(each.plannedLanduse.isPartOfRotation){
			do define_related_cells();
		}		
		aphidImmigrationIndex<-0;
		movementMask[0]<-matrix(mask1);
		movementMask[1]<-matrix(mask2);
		movementMask[2]<-matrix(mask3);
		movementMask[3]<-matrix(mask4);
		movementMask[4]<-matrix(mask5);
		movementMask[5]<-matrix(mask6);
		movementMask[6]<-matrix(mask7);
		movementMask[7]<-matrix(mask8);
										
	}
		
	reflex timestep_update { 
	//Every timestep the model does the following actions		
		ask climate{
			do update_wind;
			do update_temperature;
		}
				
	}
	
	reflex daily_update when: hour=0 {
	//Every day the model does the following actions
		ask climate{
			do update_rain;
			do update_aphid_rates;							// Update daily aphid rates (survival, repro...)
			do check_dispersal_conditions;					// Define if aphid migration is possible							
		}
		do update_landcovers;
		do updateSurf;
		do update_habitats;
		
		ask raster_landscape where(each.relatedPatch != nil and each.relatedPatch.plannedLanduse.isPartOfRotation){
			do test_if_habitat;
		}
		
		ask raster_landscape where(each.aphidPopulation>0){
			do aphid_dynamics;								// Update aphid population
			if min(myself.climate.meanTemp5LastDays) > 10.0{
				do aphid_background_mortality;				
			}
					
		}
		
		ask raster_landscape where(each.aphidPopulation>0){		
			do aphid_dispersal;								// Update aphid movements
			if localDispersers > 0 and myself.climate.hasFavourableTemperature {
				do disperse_local_flyers;			
			}

		}

		ask raster_landscape where(each.aphidPopulation>0){		
			do update_cell;
		}		
		ask raster_landscape where(each.totalArrivingWingedAphids>0.0){
			do place_winged_aphids_in_pop;
		}

		do aphid_immigration;								// Update aphid immigration		

		ask patchesList {
			do update_aphid_density();
		}
		do print_densities;
		do print_conditions;		
			
	}
	
	// VALIDATION SURFACES
	action updateSurf{
		totalWheat<-0.0;
		totalCorn<-0.0;
		totalWheatVolunteer<-0.0;
		totalSorghum<-0.0;			
		ask Landscape(0).patchesList{
			if currentLandcover.name = "CoverWheat"{
				Landscape(0).totalWheat <-Landscape(0).totalWheat + size;
			}
			if currentLandcover.name = "CoverCorn"{
				Landscape(0).totalCorn <-Landscape(0).totalCorn + size;
			}
			if currentLandcover.name = "CoverWheatVolunteer"{
				Landscape(0).totalWheatVolunteer <-Landscape(0).totalWheatVolunteer+ size;
			}
			if currentLandcover.name = "CoverSorghum"{
				Landscape(0).totalSorghum <-Landscape(0).totalSorghum+ size;	
			}
		}
		save [julianDay,year,totalWheat,totalCorn,totalWheatVolunteer,totalSorghum] to: "../models/"+"SurfacesLowPred"+string(rep)+".csv" type:"csv" rewrite:false;
	}
	
	action print_conditions{
		ask patchesList where (each.plannedLanduse.isPhenologicallyDetailed){
			save [name,size,julianDay,year,plannedLanduse,currentLandcover,phenologicalStage] to: "../models/"+"ConditionsLowPred"+string(rep)+".csv" type:"csv" rewrite:false;				
		}

	}

	action print_densities{
		ask patchesList where (each.aphidDensity >0){
			save [name,size,julianDay,year,plannedLanduse,currentLandcover,phenologicalStageIndex,aphidDensity,totalAphidDensity,minAphidDensity,maxAphidDensity,isInfested] to: "../models/"+"AphidDensitiesLowPred"+string(rep)+".csv" type:"csv" rewrite:false;	
		}
	}			
		
	action create_landuses{
	// Purpose : Create the list of Covers used by the landscape
	// In: world.coversFile 
	// Out: self.coversList
		matrix landuseMatrix <-matrix(landuseCSV);
		matrix coverStages <- matrix(cropStagesCSV);
		loop x from: 1 to: ((landuseMatrix column_at(0)) count (each!='')-1){
		// Each loop creates a cover 
			create Landuse number:1 returns: newCover{
			// Purpose : init the covers from the coversMatrix and the index of the cover in this matrix (number of the row)
				name <- string(landuseMatrix[0,x]);
				abbreviation <- string(landuseMatrix[1,x]);					
				isDynamicThroughTime <- bool(landuseMatrix[2,x]);
				if isDynamicThroughTime {
					isPartOfRotation <- bool(landuseMatrix[3,x]);
					isPhenologicallyDetailed <- bool(landuseMatrix[4,x]);
					isAnnual <- bool(landuseMatrix[5,x]);
					if isAnnual or isPartOfRotation{
						startSowingDate <- int(landuseMatrix[6,x]);
						maximumHarvestDate <- int(landuseMatrix[7,x]);	
					}
	
					if isPhenologicallyDetailed or name="WheatVolunteer" {
						baseTemperature <-int(landuseMatrix[8,x]);
						if isPartOfRotation or isAnnual {
							cropSowingDelay<-int(coverStages[1,x-1]);
							cropHarvestDelay<-int(coverStages[2,x-1]);
							int y <- 3;
							loop while: y<=(((coverStages row_at(x-1)) count (each!='')-1)){
								phenologicalStages <- phenologicalStages + coverStages[y,x-1];							
								GDDThresholds <- GDDThresholds + coverStages[y+1,x-1];
								stageResource<-stageResource + coverStages[y+2,x-1];
								phenologicalHeight<-phenologicalHeight + coverStages[y+3,x-1];	
								y<-y+4;						
							}
							nbCropStages<-length(phenologicalStages);
						}
						
					}		
					isClustered<- bool(landuseMatrix[9,x]);
				}	
				index <- x-1;										
			}
			landusesList <- landusesList + newCover;
		}
	}
		
		
	action create_landcovers{
		matrix landcoverMatrix <-matrix(landcoverCSV);
		matrix aphidCropStages <- matrix(aphidCropRates);
		loop x from: 0 to: ((landcoverMatrix column_at(0)) count (each!='')-1){
		// Each loop creates a cover 
			create Landcover number:1 returns: newCover{
			// Purpose : init the covers from the coversMatrix and the index of the cover in this matrix (number of the row)
				name <- string(landcoverMatrix[0,x]);
				abbreviation <- string(landcoverMatrix[1,x]);					
				color <- rgb([int(landcoverMatrix[2,x]),int(landcoverMatrix[3,x]),int(landcoverMatrix[4,x])]);
				baseHeight <- int(landcoverMatrix[5,x]);
				loop i from:0 to:3{
					if aphidCropStages[0,i]=name{
						loop j from:1 to:6{
							cropStageAphidRates <- cropStageAphidRates + aphidCropStages[j,i];
						} 
					}
				}					
			}							
			landcoversList <- landcoversList + newCover;
		}
	
	}
	
			
	action create_rotations{
	// Purpose : Create the list of Rotations used by the landscape
	// In: world.rotationsFile 
	// Out:self.rotationsList
		matrix rotationsMatrix<- matrix(rotationsCSV);			
		loop y from: 0 to: (rotationsMatrix column_at(0)) count (each !='')-1{
		// Each loop creates a rotation	
			create Rotation number:1 returns: newRotation{
			//Purpose : init rotations from the rotationsMatrix and the index of the rotation in this matrix (number of the row)
				list<int> phenologyTransit <- list<int>('');
				int k <- 0;
				bool startWithCover <- true;					
				name <- (string(rotationsMatrix[0,y]));
				referenceArea <- int(rotationsMatrix[1,y]);
				index <- y;
				loop x from:2 to: (rotationsMatrix row_at(y)) count (each !="")-1{
					Landuse successionCrop <- myself.landusesList first_with (each.name = (rotationsMatrix[x,y]));
					succession <- succession + successionCrop;
					phenologyDates<-phenologyDates+(successionCrop.startSowingDate-1);
					phenologyDates<-phenologyDates+successionCrop.maximumHarvestDate;
					if x = 2{
						firstDate <- successionCrop.startSowingDate;
					}
				}
				loop z from: 1 to: phenologyDates count (each !='')-1{
					if phenologyDates[z-1] > phenologyDates[z]{
						loop x from:z to: phenologyDates count (each !='')-1{
							phenologyDates[x]<-phenologyDates[x]+365;
						}
					}
				}
				duration <- int(ceil(((phenologyDates[phenologyDates count (each!='')-1]) - (firstDate)) /365));
			}			
			rotationsList <-rotationsList + newRotation;
		}
	
	}
		
	action create_patches{
	// Purpose : Create the list of Patches used by the landscape
	// In: world.landscapeFile
	// Out: self.patchesList	
		if hasInputRotations =false{				
			create Patch from: landscapeGIS with: [readLanduse::string(read("Cover")), size::int(read("Area")),iD::int(read("FID"))] returns:patchesCreated{ //
//			create Patch from: landscapeGIS with: [readLanduse::string(read("Year_2")), size::int(read("Shape_Area")),iD::int(read("OBJECTID"))] returns:patchesCreated{
				Landuse sameLanduse <- myself.landusesList first_with (each.name=self.readLanduse);  
				currentLandcover <- myself.landcoversList first_with (each.name="Cover"+sameLanduse.name);
				plannedLanduse <- sameLanduse;
				GDDCumulated<-0.0;
				isAvailableAtInit <- false;
				isInitialPatch <- false;
				hasRotationAssigned<- false;
				isVolunteer<-false;					
			}
			patchesList <-patchesCreated;
			do assign_rotations;	
			do setup_phenology;				
		}
		else {
			create Patch from: landscapeGIS with: [readLanduse::string(read("Cover")),readRotation::string(read("Rotation")),successionIndex::int(read("SucIndex"))] returns:patchesCreated{
				Landuse sameLanduse <- myself.landusesList first_with (each.name=self.readLanduse);  
				rotation <- myself.rotationsList first_with (each.name=self.readRotation);  
				currentLandcover <- myself.landcoversList first_with (each.name="Cover"+sameLanduse.name);
				plannedLanduse <- sameLanduse;
				GDDCumulated<-0.0;
				isAvailableAtInit <- false;		
			}
			patchesList<-patchesCreated;
			do setup_phenology;					
		}
	}
		
		
	action assign_rotations{
	// Purpose: assign a rotation to the agricultural patch
	// In: attribRota, coversList, rotationsList,  
	// Out: self.successionIndex, self.rotation, self.plannedCover  
	// InOut: self.actualCover
		loop i from: 0 to: (rotationsList count (each!='')-1) {	// For each rotation
			ask Rotation where (each.index=i){
				ask any(Landscape(0).patchesList where ((( each.plannedLanduse = any (self.succession))) and(each.hasRotationAssigned=false))){
					self.rotation <- myself;
					self.isInitialPatch <-true;
					self.hasRotationAssigned<-true;
					do assign_initial_landuse(rotation);				
					plannedLanduse <- rotation.succession[successionIndex];	
				}
			}
		}
		loop i from: 0 to: (rotationsList count (each!='')-2) {
			Rotation refRotation <- rotationsList[i];
			ask any(Patch where(each.rotation= refRotation)){
				refRotation.areaAssigned <- int(size);
				if rotation.succession count (each.isClustered)>0{
					loop while: (refRotation.areaAssigned < refRotation.referenceArea-(refRotation.referenceArea*2)/100){														
						Patch Neighbor <- Landscape(0).patchesList where((each.hasRotationAssigned=false) and(each.plannedLanduse.isPartOfRotation)and(each.size <((refRotation.referenceArea+(refRotation.referenceArea*5)/100))-refRotation.areaAssigned)) closest_to self;// {
						ask Neighbor{
							rotation <- myself.rotation;
							hasRotationAssigned<-true;
							do assign_initial_landuse(rotation);
							plannedLanduse <- rotation.succession[successionIndex];															
							refRotation.areaAssigned <- int(refRotation.areaAssigned + size);
						}
					}
				}
				else{
					loop while: (refRotation.areaAssigned < refRotation.referenceArea-(refRotation.referenceArea*2)/100){						
						ask any (Patch where((each.hasRotationAssigned=false) and(each.plannedLanduse.isPartOfRotation))){
							rotation <- any (Rotation where (each.index=(i)));
							hasRotationAssigned<-true;
							do assign_initial_landuse(rotation);							
							plannedLanduse <- rotation.succession[successionIndex];						
							refRotation.areaAssigned <- int(refRotation.areaAssigned + size);	
						}
					}
				}
			}
		}

		ask (Patch where((each.hasRotationAssigned=false) and(each.plannedLanduse.isPartOfRotation))){
			rotation <- any (Rotation where (each.index=(Landscape(0).rotationsList count(each!='')-1)));
			hasRotationAssigned<-true;
			rotation.areaAssigned <- int(rotation.areaAssigned + size);				
			do assign_initial_landuse(rotation);
			plannedLanduse <- rotation.succession[successionIndex];	
		}

	}	
		
	action setup_phenology{
	// Purpose: Define the phenology of each cover at the initialization
	// In : coversList, self.plannedCover, world.day
	// Out: self.actualCover
		ask Patch where (each.plannedLanduse.isPartOfRotation){
			if plannedLanduse.isAnnual {					
				currentLandcover <- myself.landcoversList first_with (each.name="Cover"+plannedLanduse.name);				// All patches are first assigned a naked landcover
				phenologicalStageIndex <- -1;										
				if isAvailableAtInit {																		// Assigns an initial landcover for fields where the crop should already be grown
					currentLandcover<- myself.landcoversList first_with (each.name="Cover"+plannedLanduse.name);				
					phenologicalStage<-'';						
				} 
				if plannedLanduse.isPhenologicallyDetailed and isAvailableAtInit =false {					// Calculates sowing and harvest delays for fields that will grow a detailed crop
					fieldSowingDelay <- rnd(plannedLanduse.cropSowingDelay);
					fieldHarvestDelay <- rnd(plannedLanduse.cropHarvestDelay);
					dayHarvest<--1;						
				}
				if plannedLanduse.isPhenologicallyDetailed and isAvailableAtInit {							// All crop availabale at init and phenologically detailed -> Flowering stage
					phenologicalStageIndex <- 3;
				}						
			}
			else{
				phenologicalStageIndex <- -1;
			}
		}			
	}			

	action update_landcovers{
	// Purpose: updates the actualCover of each patch according to it's availability. 
	//          When the cover is no longer available (actualCover = CoverBareGround), 
	//          then the coverYear changes according to the patch rotation.
	// In: rainedToday, coversList, self.plannedCover, world.day 
	// Out:  
	// InOut: self.actualCover, self.randomSowingDelay, self.successionIndex
		ask Patch where (each.plannedLanduse.isPartOfRotation){
			if (julianDay=228 and isVolunteer){
				currentLandcover <- myself.landcoversList first_with (each.name="CoverWheatVolunteer");
				isVolunteer <- false;
			}
			if (julianDay=344 and currentLandcover.name="CoverWheatVolunteer"){
				currentLandcover <- myself.landcoversList first_with (each.name="CoverBareGround");
			}										
			if (plannedLanduse.isDynamicThroughTime) = true{
				if plannedLanduse.isPartOfRotation = true{
					if plannedLanduse.isPhenologicallyDetailed = true and isAvailableAtInit=false{
						//////// CALCULATE DELAYS ///////////						
						int randomSowing<-plannedLanduse.startSowingDate+fieldSowingDelay;				
						if randomSowing >365{
							randomSowing<-randomSowing-365;
						}
						/////////  SOW IF NOT RAINING //////////
						if julianDay=randomSowing and myself.climate.isRainingToday = false{
							currentLandcover <-  myself.landcoversList first_with (each.name="CoverBareGround");
							GDDCumulated<-0.0;
							phenologicalStageIndex<-0;
							phenologicalStage <- "Sowed";
							isVolunteer<-false;				
						}
						////// IF RAINING DELAY TO TOMORROW /////						
						if julianDay=randomSowing and  myself.climate.isRainingToday = true {
							fieldSowingDelay <- fieldSowingDelay +1;
						}
						////// IF CROP IS AT ANY STAGE /////////						
						if (phenologicalStageIndex >= 0) {	
							/// CALCULATE GDD //
							if (((myself.climate.minDayTemp+myself.climate.maxDayTemp)/2)-plannedLanduse.baseTemperature) > 0{
								GDDCumulated <- GDDCumulated +  (((myself.climate.minDayTemp+myself.climate.maxDayTemp)/2)-plannedLanduse.baseTemperature) ;
							}							
							/// UPDATE STAGES ///
							if (phenologicalStageIndex<(plannedLanduse.nbCropStages)) and (GDDCumulated >= plannedLanduse.GDDThresholds[phenologicalStageIndex]){
								currentLandcover<- myself.landcoversList first_with (each.name="Cover"+plannedLanduse.name);					
								phenologicalStage <- plannedLanduse.phenologicalStages[phenologicalStageIndex];
								//DEFINE HARVEST DATE//
								if phenologicalStage = "Harvestable"{
									dayHarvest <- julianDay + fieldHarvestDelay;
									if dayHarvest > 365{
										dayHarvest <- dayHarvest - 365;
									}
								}
								phenologicalStageIndex <- phenologicalStageIndex +1;
							}				
						}
						// HARVEST //
						if (julianDay=dayHarvest)and myself.climate.isRainingToday {
							dayHarvest <- dayHarvest +1;
						}
						
						if (julianDay=dayHarvest) and myself.climate.isRainingToday=false {
							if currentLandcover.name='CoverWheat'{
								currentLandcover <- myself.landcoversList first_with (each.name="CoverBareGround");
								isVolunteer<-true;
							}
							else{
							currentLandcover <-  myself.landcoversList first_with (each.name="CoverBareGround");								
							}
							phenologicalStage <-"";
							phenologicalStageIndex<--1;								
							successionIndex <- successionIndex +1;
							if successionIndex > (length(rotation.succession))-1{
								successionIndex<-0;
							}

							plannedLanduse <- rotation.succession[successionIndex];						
							if plannedLanduse.isPhenologicallyDetailed {
								fieldSowingDelay <- rnd(plannedLanduse.cropSowingDelay);
								fieldHarvestDelay <- rnd(plannedLanduse.cropHarvestDelay);
								dayHarvest<--1;									
							}								
				
						}						
						if (julianDay=plannedLanduse.maximumHarvestDate and currentLandcover.name!="CoverBareGround"){
							if currentLandcover.name='CoverWheat'{
								currentLandcover <- myself.landcoversList first_with (each.name="CoverBareGround");
								isVolunteer <- true;
							}
							else{
							currentLandcover <-  myself.landcoversList first_with (each.name="CoverBareGround");								
							}							
							phenologicalStage <-"";
							phenologicalStageIndex<--1;								
							successionIndex <- successionIndex +1;
							if successionIndex > (length(rotation.succession))-1{
								successionIndex<-0;
							}

							plannedLanduse <- rotation.succession[successionIndex];									
							if plannedLanduse.isPhenologicallyDetailed {
								fieldSowingDelay <- rnd(plannedLanduse.cropSowingDelay);
								fieldHarvestDelay <- rnd(plannedLanduse.cropHarvestDelay);
								dayHarvest<--1;									
							}								
												
						}							
					}
					else{
						if julianDay=plannedLanduse.startSowingDate{
							currentLandcover<- myself.landcoversList first_with (each.name="Cover"+plannedLanduse.name);
							phenologicalStage <-"";
							phenologicalStageIndex <- -1;
							isVolunteer<-false;							
						}
						if julianDay= plannedLanduse.maximumHarvestDate and currentLandcover.name="Cover"+plannedLanduse.name{
							currentLandcover <-  myself.landcoversList first_with (each.name="CoverBareGround");						
							phenologicalStage <-"";
							phenologicalStageIndex <- -1;								
							successionIndex <- successionIndex + 1;
							if successionIndex > (length(rotation.succession))-1{
								successionIndex<-0;
							}
							isAvailableAtInit<- false;
							plannedLanduse <- rotation.succession[successionIndex];									
							if plannedLanduse.isPhenologicallyDetailed {
								fieldSowingDelay <- rnd(plannedLanduse.cropSowingDelay);
								fieldHarvestDelay <- rnd(plannedLanduse.cropHarvestDelay);
								dayHarvest<--1;																									
							}																			
						}	
					}
				}		
			}	
		}
		ask patchesList{
			if currentLandcover.name = "Cover"+plannedLanduse and plannedLanduse.isPhenologicallyDetailed and isAvailableAtInit=false{
				patchCoverHeight <- plannedLanduse.phenologicalHeight[phenologicalStageIndex-1];
			}
			else {
				patchCoverHeight <- currentLandcover.baseHeight;
			}
		}							
	}
		
	action update_habitats{
			availableHabitats <- raster_landscape  where(each.habitat=1);
			save [julianDay,year,length(availableHabitats)] to: "../models/"+"HabitatsRef1"+".csv" type:"csv" rewrite:false; 							
	}	

	action create_aphid_test{
		raster_landscape testAphids <- any(raster_landscape where (each.relatedPatch != nil and ((each.relatedPatch.plannedLanduse.isPhenologicallyDetailed and each.relatedPatch.isAvailableAtInit)or(each.relatedPatch.currentLandcover.name='CoverWheatVolunteer'))));	
		ask testAphids {
			numberAdults  <- 10;
			aphidPopulation <- 10;
		}
	}
		
	action aphid_immigration{ 
	// Calculate the number of immigrating alates and place them in the landscape depending on suction trap data	
		int dailyCapturedFlyers <- int(dailyImmigrantData[0,aphidImmigrationIndex]);			// Number of captured alates in the trap for the current day
		int totalFieldArea <- 0;
		ask Landscape(0).patchesList where ((each.plannedLanduse.isPhenologicallyDetailed  and each.phenologicalStageIndex > 0)or(each.currentLandcover.name='CoverWheatVolunteer')){
			if currentLandcover.name="CoverWheatVolunteer" or plannedLanduse.stageResource[phenologicalStageIndex-1] {
				totalFieldArea <- int(totalFieldArea + size);							// Here calculate the total area covered by potential habitats
			}
		}
		totalFieldArea <- int(totalFieldArea/10000);
		if totalFieldArea >0{
			int numberImmigrants <- totalFieldArea*237*dailyCapturedFlyers;					// Estimate the number of immigrants with the Morgan (2000) relation
			if numberImmigrants >0{
				if length(availableHabitats)>0{
					list<raster_landscape> edgeCells <- availableHabitats where (each.isEdgeCell);
					list<raster_landscape> centerCells <- availableHabitats where (each.isEdgeCell=false);			
					list<float> cellProbs;
					loop i from:0 to: length(edgeCells)-1{
						cellProbs<-cellProbs + (10/((10*length(edgeCells)+length(centerCells))));
					}
					loop j from:0 to: length(centerCells)-1{
						cellProbs<-cellProbs + (1/((10*length(edgeCells)+length(centerCells))));				
					}
					list<float> totalCellProbs;
					loop i from:0 to: length(cellProbs)-1{
						float k <- 0.0;
						loop j from:0 to:i{
							k<- k + cellProbs[j];
						}
						if i=length(cellProbs)-1{
							k<-1.0;
						}
						totalCellProbs<-totalCellProbs + k;
					}
					list<raster_landscape> trueHabitats <- edgeCells + centerCells;
					int i <- length(trueHabitats)-1;
					loop while: i>=0{
						int cellArrivants <-0;
						if i=0{
							cellArrivants <- numberImmigrants;
						}
						else{
							cellArrivants <- binomial(numberImmigrants,cellProbs[i]/totalCellProbs[i]);
						}
			 			if cellArrivants >0{
			 				ask trueHabitats[i]{
			 					numberAdults  <- numberAdults +cellArrivants;
								aphidPopulation <- aphidPopulation +cellArrivants;
			 				}
			 				numberImmigrants <- numberImmigrants - cellArrivants;			 								
						}
						if numberImmigrants=0{
							i<--1;
						}
						else{
			 				i <- i-1;							
						}
			 		}	 			
				}			
			}				
		
		
		}
		aphidImmigrationIndex <- aphidImmigrationIndex+1;
		int	aphidImmigrationIndexMax <-((dailyImmigrantData column_at(0)) count (each !='')-1);
		if aphidImmigrationIndex > aphidImmigrationIndexMax {aphidImmigrationIndex <-1;}		
		
	}


		
	action calculate_movement_probabilities{		// Calculates each movement circle threshold, then the average circle for moving the aphids within a group and finally the 
													// probabilities of being in each group according to the defined gamma distribution describing aphid dispersal		
		loop i from: 0 to: (number_movement_circles-2){
			distanceGroupProbabilities[i] <- pgamma((grid_cell_size/2)+(i*grid_cell_size),1.75513,0.03487638);			
		}	
		loop k from:0 to:(number_movement_circles-2){
			if k=0{
				aphidMovementDistProbs[k]<-distanceGroupProbabilities[k];
			}
			else{
				aphidMovementDistProbs[k]<-distanceGroupProbabilities[k]-distanceGroupProbabilities[k-1];				
			}			
		}
		distanceGroupProbabilities<-distanceGroupProbabilities+1;
		aphidMovementDistProbs<-aphidMovementDistProbs +(1-sum(aphidMovementDistProbs));		
	}				

}
		
	
species Climate{
	float windSpeed; 	//Speed of the wind
	float windDirection;// Direction of the wind (360° degrees)
	int tempIndex; 		//Index in the temperature matrix
	int rainIndex;		// Index in the rain matrix
	bool isRainingToday;	//true if it has rained today ?
	float landscapeTemperature; // Temperature in the landscape (as data from the weather station)
	matrix<float> dayTemperatures <- 0.0 as_matrix{24/timestepHours,1} ; // List of the temperatures of the day (every timestep)
	float meanDayTemp; 		//mean temperature of the day
	float minDayTemp;		//min temperature of the day
	float maxDayTemp;		//max temperature of the day
	matrix dailyRainData <- 0.0 as_matrix{2,1826}; 	//input csv as a matrix
	matrix windTempData;	//input csv as a matrix 
	list<float> meanTemp5LastDays <-[0.0,0.0,0.0,0.0,0.0];

	float devApNymphs;		//daily development rate of apterous nymphs
	float devAlNymphs;		//daily development rate of alate nymphs
	float devAdults;		//daily development rate of adults
	float aphidReproductionRate;	//daily reproduction rate
	float aphidLethalDailyDegrees;	//daily degrees below 2.8°c increasing mortality
	float aphidSurvivalRate;		//daily survival rate
	bool hasFavourableMigration;	//Weather conditions are favourable for migration flight 
	bool hasFavourableTemperature;			//
	bool hasWindInfluence;
	float meanDayWind;				// Mean wind speed value of the day
	matrix<float> dayWindSpeeds <- 0.0 as_matrix{24/timestepHours,1}; // List of all daily windspeed values.		
				
	init{
	// Purpose: initialize the climate attributes
	// Out: self
		dailyRainData <- matrix(rainCSV); //Matrix of the rainFile
		windTempData <- matrix(climateCsv); //Matrix of the climatFile (wind+temperatures)
		isRainingToday <- false; // it has not rained today
		rainIndex<-1; //Index in the rain matrix								
		tempIndex <- 1; //Index in the climat matrix
		
	}
		
	action update_wind {
	// Purpose: Update the wind speed and direction at every step
	// In: self.windTempData
	// Out: self.windSpeed, self.windDirection			
		windSpeed <- float(windTempData[5,tempIndex]);
		windDirection <- float(windTempData[4,tempIndex]);
		int i<-int(hour/timestepHours);
		dayWindSpeeds[i,0]<-windSpeed;
		if hour = 0{				
			meanDayWind <- mean(dayWindSpeeds);
		}	
	}
		
	action update_temperature {
	// Purpose: Update the temperature at every step	
	// In: self.windTempData
	// Out: self.landscapeTemperature, self.meanDayTemp	
	// InOut: self.tempIndex																							
		landscapeTemperature <- float(windTempData[2,tempIndex]);
		int i<- int(hour/timestepHours);
		dayTemperatures[i,0] <-(landscapeTemperature);
		if hour=0{
			if hasHourlyWeatherData{
				meanDayTemp <- mean(dayTemperatures);
				minDayTemp <- min(dayTemperatures);
				maxDayTemp<- max(dayTemperatures);			
			}				
			else{
				minDayTemp <- dayTemperatures[0,0];
				maxDayTemp<-dayTemperatures[1,0];
				meanDayTemp <- dayTemperatures[2,0];										
			}
		}				
		tempIndex <- tempIndex + 1;
		meanTemp5LastDays[4]<-meanTemp5LastDays[3];
		meanTemp5LastDays[3]<-meanTemp5LastDays[2];
		meanTemp5LastDays[2]<-meanTemp5LastDays[1];
		meanTemp5LastDays[1]<-meanTemp5LastDays[0];			
		meanTemp5LastDays[0]<-meanDayTemp;		
		int tempIndMax <- ((windTempData column_at(0)) count (each !='')-1);			
		// Number of lines in this matrix (ends september 14th in order to be able to loop)	
		if tempIndex > tempIndMax {tempIndex <- 1;}
	}
		
	action update_rain{
	// Purpose: Update the rain every day
	// In: self.dailyRainData
	// Out: self.rainedToday
	// InOut: self.rainInd	
		int raining <- int(dailyRainData[1,rainIndex]);
		if raining = 1{isRainingToday <- true;} else {isRainingToday<-false;}
		rainIndex <- rainIndex +1;
		int	rainIndMax <-((dailyRainData column_at(0)) count (each !='')-1);
		if rainIndex > rainIndMax {rainIndex <-1;}
	}
		
	action update_aphid_rates{
	//Purpose : Update the different temperature/weather dependant rates affecting aphid dynamics
			
////////////////////////			// Development Rates //					/////////////////////////////////////////////////////////////
		aphidLethalDailyDegrees <- 0.0;
		float apNyDevMin <- -0.015 + (0.291 /(1+exp(-0.138*(minDayTemp-16.911)))); 	//update apNymphs rate @min daily temperature
		float apNyDevMean <- -0.015 + (0.291 /(1+exp(-0.138*(meanDayTemp-16.911))));	//update apNymphs rate @mean daily temperature
		float apNyDevMax <- -0.015 + (0.291 /(1+exp(-0.138*(maxDayTemp-16.911))));		//update apNymphs rate @max daily temperature
		float alNyDevMin <- apNyDevMin /1.5;										//update alNymphs rate @min daily temperature
		float alNyDevMean <- apNyDevMean /1.5;										//update alNymphs rate @mean daily temperature
		float alNyDevMax <- apNyDevMax /1.5;										//update alNymphs rate @max daily temperature
		float adDevMin <- 0.0193 + 0.0039 * minDayTemp;								//update adult rate @min daily temperature
		float adDevMean <- 0.0193 + 0.0039 * meanDayTemp;							//update adult rate @mean daily temperature
		float adDevMax <- 0.0193 + 0.0039 * maxDayTemp;								//update adult rate @max daily temperature			
		//Extreme temperatures lead to a development rate of 0
		if (minDayTemp <=5.8)or (minDayTemp>=30){
			apNyDevMin<-0.0;
			alNyDevMin<-0.0;
			adDevMin<-0.0;
		}
		if (meanDayTemp <=5.8)or (meanDayTemp>=30){
			apNyDevMean<-0.0;
			alNyDevMean<-0.0;
			adDevMean<-0.0;
		}
		if (maxDayTemp <=5.8)or (maxDayTemp>=30){
			apNyDevMax<-0.0;
			alNyDevMax<-0.0;
			adDevMax<-0.0;
		}						
			
		devAlNymphs <- (alNyDevMin+alNyDevMean+alNyDevMax)/3;		//daily alNymph dev rate using simpsons rule
		devApNymphs <- (apNyDevMin+apNyDevMean+apNyDevMax)/3;		//daily apNymph dev rate using simpsons rule			
		devAdults <- (adDevMin+adDevMean+adDevMax)/3;		 		//daily adult dev rate using simpsons rule
			
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////			
			
			
//////////////////////////////////			// Reproduction rates //              ////////////////////////////////////////////////
		float reproductionRateMin <- -0.036 + (5.825 /(1+exp(-0.319*(minDayTemp-12.03))));   //daily repro @min daily temperature
		float reproductionRateMean <- -0.036 + (5.825 /(1+exp(-0.319*(meanDayTemp-12.03)))); //daily repro @mean daily temperature
		float reproductionRateMax <- -0.036 + (5.825 /(1+exp(-0.319*(maxDayTemp-12.03))));   //daily repro @max daily temperature
		
		//Extreme temperatures
		if (minDayTemp <=5.8)or (minDayTemp>=30){
			reproductionRateMin<-0.0;
		}
		if (meanDayTemp <=5.8)or (meanDayTemp>=30){
			reproductionRateMean<-0.0;
		}
		if (maxDayTemp <=5.8)or (maxDayTemp>=30){
			reproductionRateMax<-0.0;
		}			

		aphidReproductionRate <-(reproductionRateMin+ reproductionRateMin + reproductionRateMax)/3; //Simpsons rule
		
		
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////		
		
		
////////////////////////////////////			// Survival Rate //          //////////////////////////////////////////////////////////
		
		if ((minDayTemp+maxDayTemp)/2) <= 2.8{ //calculate degree days below 2.8
			aphidLethalDailyDegrees <-2.8 - ((minDayTemp+maxDayTemp)/2);
		}
			
		aphidSurvivalRate <- 0.9511 - 0.0173  * aphidLethalDailyDegrees; //calculate the survival rate
			
		// Extreme temperatures
		if (((minDayTemp+meanDayTemp+maxDayTemp)/3) >= 30)or(((minDayTemp+meanDayTemp+maxDayTemp)/3) <= -4){ //if temperature >=30°c or <=-4°C, mortality rate is of 0.5
			aphidSurvivalRate <- 0.5;
		}
		
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////			
	}
		
	action check_dispersal_conditions{
		hasFavourableMigration <- false;			
		int seasonLowerThreshold<-9;			
		if (julianDay >=79) and(julianDay<=172){
			seasonLowerThreshold <- 16;
		}
		if (julianDay >=173)and (julianDay<=265){
			seasonLowerThreshold <- 13;
		}
		if (meanDayWind<11.0) and (((minDayTemp+meanDayTemp+maxDayTemp)/3) > seasonLowerThreshold){
			hasFavourableMigration <- true;
		}
			
			////////////////////////////: TEST SORTIE CONDITIONS///////////////////////			
		if (meanDayWind<11.0){
			hasWindInfluence <-true;
		}
		else {
			hasWindInfluence <-false;
		}
		if (((minDayTemp+meanDayTemp+maxDayTemp)/3) > seasonLowerThreshold){
			hasFavourableTemperature<-true;
		}
		else{
			hasFavourableTemperature<-false;
		}
	/////////////////////////////////////////////////			
	}	
}
	

	
species Patch{
	string readRotation; 	// Rotation read in shapefile if info already in the input
	string readLanduse;		// Cover read in the shapefile
	int iD;
	float size; // Patch size in m²		
	bool isInitialPatch;
	bool hasRotationAssigned;
	Rotation rotation; 		// Rotation assigned to the patch
	bool isAvailableAtInit; // When we init, we have no info on DD so we use the dates for harvesting
	int successionIndex;	// Index of the cover in the succession of the rotation	
	Landuse plannedLanduse; 	// Cover that will be grown this year
	Landcover currentLandcover; 		// Cover that is actually on the field		
	int fieldSowingDelay; 	// Delay due to human intervention for sowing
	int fieldHarvestDelay; // Delay due to human intervention for harvest		
	int dayHarvest;		
	float GDDCumulated; // degree days absorbed by the plant
	int phenologicalStageIndex;
	string phenologicalStage;		  
	int patchCoverHeight;
	int aphidDensity;
	int totalAphidDensity;
	int minAphidDensity;
	int maxAphidDensity;
	bool isInfested;
	list<raster_landscape> relatedCells;
	
	bool isVolunteer;
	
	action define_related_cells{
		relatedCells <- raster_landscape where (each.relatedPatch = self) ;		
	}
													
	action assign_initial_landuse (Rotation assignedRotation){
		int k <- rnd(assignedRotation.duration);
		int successionPosition <- 0;
		if assignedRotation.firstDate > julianDay{
			k<-k+1;
		}
		k<-julianDay+k*365;
		loop x from:0 to: assignedRotation.phenologyDates count (each !="")-1{
			if assignedRotation.phenologyDates[x] < k{
				successionPosition <- successionPosition +1;
			}
		}
		successionIndex<-int(successionPosition/2);
		if successionPosition = (assignedRotation.phenologyDates count (each !="")){
			successionIndex<-0;
		}
		if successionPosition mod 2 = 0{
			isAvailableAtInit<-false;
		}
		else{
			isAvailableAtInit<-true;
		}	
	}	

	action update_aphid_density{
		totalAphidDensity<-0;		
		aphidDensity <- 0;
		list<int> aphidDensities;
		if length(relatedCells) >0{
			loop i over: relatedCells where (each.aphidPopulation>0){
				aphidDensities<-aphidDensities+i.aphidPopulation;
				self.totalAphidDensity <- self.totalAphidDensity + i.aphidPopulation;				
			}
			aphidDensity <- int(totalAphidDensity/length(relatedCells));   // Aphid density per 900 m²
			minAphidDensity<-min(aphidDensities);
			maxAphidDensity<-max(aphidDensities);									
		}
		if aphidDensity >= 2250000{
			isInfested <-true;
		}
		else{
			isInfested <- false;
		}		
	}
					
	aspect base {draw shape  color: currentLandcover.color depth: patchCoverHeight;}
}
	
species Rotation{
	int duration; // Number of years covered by the rotation
	list<Landuse> succession; //The succession of covers in the rotation
	int referenceArea; //Area that has to be assigned to the rotation
	int index; //Index of the rotation in world.Landscape(0).rotationsList
	           // (used when creating the attribution matrix)
	int firstDate;           
	int areaAssigned;  
	list<int> phenologyDates;
	list<Landuse> initialLanduses;
	list<bool> isAvailableAtInit;         
}
	
species Landcover{
	rgb color;           // Color on the map
	string abbreviation;  // 3 Letter abreviation
	int baseHeight;          // Height of the cover in a 3D map
	bool resource; // Is the crop a resource in our population dynamics context ?
	int index;           // Index of the cover in Landscape.coversList 
	                     // (used when creating the attribution matrix)	
	list<float> cropStageAphidRates;		                     	                                                          			   
}	
	
species Landuse{
	string abbreviation;  // 3 Letter abreviation
	int index;           // Index of the cover in Landscape.coversList 
	                     // (used when creating the attribution matrix)  
	bool isDynamicThroughTime; // Static or dynamic
	bool isPartOfRotation ; // Will this landuse vary according to a rotation assigned to the patch
	bool isPhenologicallyDetailed; // Shoud a detailed phenology be considered for this landuse	                                    
	bool isAnnual; // Annual or perennial species		                                    
	bool isClustered; //Is the crop clustered in the landscape	
	int startSowingDate; //Date when it is possible to sow
	int maximumHarvestDate; //Maximum date to harvest
	int cropHarvestDelay; 
	int cropSowingDelay;
	list<string> phenologicalStages;
	list<int> GDDThresholds;
	list<bool> stageResource; 
	list<int> phenologicalHeight;
	int nbCropStages; 
	int baseTemperature; // base growth temperature of the crop										
}

grid raster_landscape file:(rasterFile){
	int aphidPopulation; // Total population of the cell
	int unfavourableHabitat;
	Patch relatedPatch; // Vector patch related to the cell
	int numberApNymphs; // Number of apterous nymphs
	int numberAlNymphs; // Number of alate nymphs
	int numberAdults; // Number of adults
	list<int> wingedAphids; // Number of produced winged aphids (0 to 4 days old)
	int nymphsIntoWinged;	// Number of nymphs becoming adults for the day			
	float cropStageEffect;	// Effect of the current crop and its stage on aphid reproduction
	int patchId;			// ID of the related patch used in arcGIS to associate cells to patch
	bool isEdgeCell;			// Defines the cells at the edge of a patch, which receive potentially more migrants
	int habitat;
	int localDispersers;
	matrix<int> potentialArrivingCells <- 0 as_matrix{15,15};
	float totalArrivingWingedAphids;

	
	geometry ring -> {circle(250) - circle(200)};		// Geometry used in aphid dispersal									
		
	init{
		numberApNymphs <-0;
		numberAlNymphs<-0;
		numberAdults <-0;
		wingedAphids<-[0,0,0,0];
		patchId <-int(grid_value);
		isEdgeCell<-false;
		habitat<-0;
		totalArrivingWingedAphids<-0.0;	
	}	
		
	action define_reference_patch{
	// associates cells to patches, in order to update crops etc...		
		relatedPatch <- any(Patch where (each.name='Patch'+patchId));  
	}
	
	action test_if_habitat{
		if (relatedPatch.plannedLanduse.isPhenologicallyDetailed  and relatedPatch.phenologicalStageIndex > 0){
			if (relatedPatch.plannedLanduse.stageResource[relatedPatch.phenologicalStageIndex-1])='true'{
				habitat<-1;
			}

			else{
				habitat<-0;			
			}
		}
		else{
			habitat<-0;			
		}
		if relatedPatch.currentLandcover.name="CoverWheatVolunteer"{
			habitat <- 1;
		}
		
	}
	
	action define_edge_cells{
	// defines cells at the border of patches, used in aphid dispersal
		list<raster_landscape> neighbours<- raster_landscape at_distance 1;
		raster_landscape differentPatch <- any(neighbours where(each.relatedPatch != self.relatedPatch));
		if differentPatch!=nil{
			if differentPatch.relatedPatch=nil{
				isEdgeCell<-false;		
			}
			else{
				isEdgeCell<-true;
			}						
		}
	}
	
	
	action aphid_dynamics{
		if aphidPopulation >0{
		// Applies to colonies already settled	 
///////////// Evaluate the status of the field in terms of habitat //////////////////////////////////
			if relatedPatch.currentLandcover.name='CoverWheatVolunteer' or (relatedPatch.phenologicalStageIndex>0 and relatedPatch.plannedLanduse.stageResource[relatedPatch.phenologicalStageIndex-1]){
				unfavourableHabitat <- 0;
			}
			else {
				unfavourableHabitat <- unfavourableHabitat + 1;
			}
			if relatedPatch.phenologicalStageIndex<=0{
				unfavourableHabitat<-unfavourableHabitat +1;
			}
				
		
////////////////				// Development //         //////////////////////////////////////////////////////
				
			if numberAdults >0 and Landscape(0).climate.devAdults> 0 {
				int agingAdults <- binomial(numberAdults,Landscape(0).climate.devAdults);
				numberAdults <- numberAdults - agingAdults; //Adults die of age					
			}				
			if numberApNymphs >0 and Landscape(0).climate.devApNymphs>0{
				int agingApNymphs <- binomial(numberApNymphs,Landscape(0).climate.devApNymphs);
				numberAdults <- numberAdults + agingApNymphs; // Ap Nymphs become adults
				numberApNymphs <- numberApNymphs - agingApNymphs; // Ap Nymphs that stay nymphs										
			}
			if numberAlNymphs >0 and Landscape(0).climate.devAlNymphs>0 {
				nymphsIntoWinged <- binomial(numberAlNymphs,Landscape(0).climate.devAlNymphs); // Al Nymphs becoming winged aphids
				numberAlNymphs<-numberAlNymphs - nymphsIntoWinged;						
			}
		
/////////////////////////////////// Proportion of winged aphids ///////////////////////////////////////////////
			float wingedProportion <- (0.002 + 0.991)/(1+exp(-0.076*(((numberAdults+numberApNymphs+numberAlNymphs)/35000)-67.416))); // Density dependant production					
			if relatedPatch.currentLandcover.name = "CoverWheatVolunteer"{
				wingedProportion <- (0.002 + 0.991)/(1+exp(-0.076*(((numberAdults+numberApNymphs+numberAlNymphs)/13000)-67.416)));
			}
			if unfavourableHabitat >0{ // If crop becomes mature (unfavourable) 100% production of winged aphids
				wingedProportion <-1.0;
			}

/////////////////////////// Effect of crop stage on reproduction ////////////////////////////////////////////////
			if unfavourableHabitat = 0{
				if relatedPatch.currentLandcover.name='CoverWheatVolunteer'{
					cropStageEffect <- 1.09;					
				}
				else{
					if relatedPatch.currentLandcover.name='CoverBareGround'{
						cropStageEffect <-1.0;
					}
					if relatedPatch.plannedLanduse.isPhenologicallyDetailed{
						cropStageEffect <- relatedPatch.currentLandcover.cropStageAphidRates[relatedPatch.phenologicalStageIndex-1];					
					}					
				}

			}
			else{
				cropStageEffect <- 1.0;
			}

////////////////////// Reproduction of aphids ///////////////////////////////////////////				
			int actualReproduction <-int(Landscape(0).climate.aphidReproductionRate*cropStageEffect);

			if numberAdults >0 and Landscape(0).climate.aphidReproductionRate>0{
				int newNymphs <- numberAdults * actualReproduction + binomial(numberAdults,(Landscape(0).climate.aphidReproductionRate*cropStageEffect)-actualReproduction);
				if newNymphs >0{
					int newAlNymphs <- 0;
					if wingedProportion = 1{
						newAlNymphs <- newNymphs;
					}
					else{
						newAlNymphs <- binomial(newNymphs,wingedProportion);							
					}
					numberAlNymphs <- numberAlNymphs + newAlNymphs; // calculate the number of Al nymphs born
					numberApNymphs <- numberApNymphs + (newNymphs-newAlNymphs); // Apterous							
				}					
			}
			
/////////////////////// Death	///////////////////////////////////////////	
			if numberAlNymphs >0{
				numberAlNymphs <- binomial(numberAlNymphs,Landscape(0).climate.aphidSurvivalRate);	
			}
			if numberApNymphs >0{
				numberApNymphs <- binomial(numberApNymphs,Landscape(0).climate.aphidSurvivalRate);
			}
			if numberAdults >0{
				numberAdults <- binomial(numberAdults,Landscape(0).climate.aphidSurvivalRate);	
			}	
										
			if unfavourableHabitat > 0{
				numberAlNymphs <- int(numberAlNymphs * ((3-unfavourableHabitat)/3));
				numberApNymphs <- int(numberApNymphs * ((3-unfavourableHabitat)/3));
				numberAdults <- int(numberAdults * ((3-unfavourableHabitat)/3));																
			}
						
///////////////////// Visual representation of the colony ////////////////////////////////////						
			aphidPopulation <- numberAlNymphs+numberApNymphs+numberAdults; // 
			if aphidPopulation >0.0{
				color <- rgb('red');
			}
			else{
				color<- rgb('transparent');
			}				
		}		
	}
		
	action aphid_dispersal{
		if aphidPopulation > 0{
			
////////////// If migration possible, all winged aphids fly out of the model //////////////////////////////			
			if Landscape(0).climate.hasFavourableMigration = true{					
				wingedAphids <- [0,0,0,0];										
			}
			
			
///////////////  If not, local dispersion according to movement distances /////////////////////////////	
			else {
				localDispersers <- round(wingedAphids[3]);
				wingedAphids[3] <- wingedAphids[2];
				wingedAphids[2] <- wingedAphids[1];
				wingedAphids[1] <- wingedAphids[0];
				wingedAphids[0] <- nymphsIntoWinged;					
			}	
		}
	
	}
	
	action disperse_local_flyers{
		list<int> circle_dispersers<-list_with(Landscape(0).number_movement_circles,0.0); 
		matrix<float> cellMovementMatrix <- 0  as_matrix {Landscape(0).numberMovementCells,Landscape(0).numberMovementCells};	

				
		int i<-length(circle_dispersers)-1;
		loop while: (i>=0){
			if i=0{
				circle_dispersers[i]<-localDispersers;
			}
			else {
				circle_dispersers[i]<-binomial(localDispersers,(Landscape(0).aphidMovementDistProbs[i]/Landscape(0).distanceGroupProbabilities[i]));
				localDispersers<-localDispersers-circle_dispersers[i];				
			}
			i <- i -1;
			if localDispersers = 0{
				i <- -1;
			}			
		}
		
		
		loop k from:0 to: Landscape(0).numberMovementCells-1{
			loop g from:0 to: Landscape(0).numberMovementCells-1{
				list<int> cellCoordinates <- action_get_absolute_coordinates(self.grid_x,self.grid_y,k,g, Landscape(0).landscape_grid_size_X, Landscape(0).landscape_grid_size_Y, Landscape(0).numberMovementCells);
				raster_landscape targetCell <- raster_landscape[cellCoordinates[0],cellCoordinates[1]];
				 self.potentialArrivingCells[k,g]<-targetCell.habitat;
			}
		}	


		loop j from:0 to: length(circle_dispersers)-1{
			if circle_dispersers[j]>0{
				matrix<float> circle_cells<-matrix<int>(Landscape(0).movementMask[j])*potentialArrivingCells;
				int number_available_cells <- sum(circle_cells);
				if number_available_cells = 0{
					if j > 0{
						circle_dispersers[j-1]<-circle_dispersers[j-1]+circle_dispersers[j];
					}
				}
				else{
					float indiv_per_cell <- circle_dispersers[j]/number_available_cells;
					circle_cells<-circle_cells*indiv_per_cell;
					cellMovementMatrix<-cellMovementMatrix+circle_cells;					
				}
			}
		}
		
		loop k from:0 to: Landscape(0).numberMovementCells-1{
			loop g from:0 to: Landscape(0).numberMovementCells-1{
				list<int> cellCoordinates <- action_get_absolute_coordinates(self.grid_x,self.grid_y,k,g, Landscape(0).landscape_grid_size_X, Landscape(0).landscape_grid_size_Y, Landscape(0).numberMovementCells);
				raster_landscape targetCell <- raster_landscape[cellCoordinates[0],cellCoordinates[1]];
				 targetCell.totalArrivingWingedAphids<-targetCell.totalArrivingWingedAphids+cellMovementMatrix[k,g];
			}
		}			

	}
	
	action place_winged_aphids_in_pop{
		numberAdults<-numberAdults+round(totalArrivingWingedAphids);
		totalArrivingWingedAphids<-0.0;		
	}
	
	action aphid_background_mortality{
		numberApNymphs<-int(numberApNymphs*0.7);
		numberAlNymphs<-int(numberAlNymphs*0.7);
		numberAdults<- int(numberAdults*0.7);		
	}	
	
	action update_cell{
		if aphidPopulation >0.0{
			color <- rgb('red');
		}
		else{
			color<- rgb('transparent');
		}			
	}
	
	list<int> action_get_absolute_coordinates(int cellX,int cellY,int matrixX,int matrixY, int landscape_grid_size_X, int landscape_grid_size_Y, int window_size){
		int xAbs <- cellX+matrixX-((window_size-1)/2);
		int yAbs <- cellY+matrixY-((window_size-1)/2);
		if (xAbs <0){
			xAbs<-landscape_grid_size_X+xAbs;
		}
		if (xAbs> (landscape_grid_size_X-1)){
			xAbs<-xAbs-landscape_grid_size_X;
		}
		if (yAbs <0){
			yAbs<-landscape_grid_size_Y+yAbs;
		}
		if (yAbs> (landscape_grid_size_Y-1)){
			yAbs<-yAbs-landscape_grid_size_Y;
		}
		return [xAbs,yAbs];
				
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

