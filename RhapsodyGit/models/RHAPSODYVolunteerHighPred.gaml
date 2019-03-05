/**
 *  CADDYmodel
 *  Author: Dynafor
 *  Description: 
 */

model CADDYmodel

import 'RHAPSODY_LauncherHighPred.gaml'
import 'ATLAS_RHAPSODYVolunteerHighPred.gaml'


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
		if differentPatch.relatedPatch=nil{
			isEdgeCell<-false;		
		}
		else{
			isEdgeCell<-true;
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
		numberApNymphs<-int(numberApNymphs*0.5);
		numberAlNymphs<-int(numberAlNymphs*0.5);
		numberAdults<- int(numberAdults*0.5);		
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
	

