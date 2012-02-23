float repConst = 1.0;
float migConst = 1/10.0;
float fluConst = 0.0;
float adhConst = 0.1;

float radius = 10.0;
float k = 1.0;
float epsilon = 0.01;
float pi = 3.14159;
ArrayList cells;
int numCells = 20;
int numMigrants = 20;
PVector nurseLocation = new PVector(600,400);
PVector borderLocation = new PVector(200,400);
float nurseDisperse = 100.0;
float borderDisperse = 10.0;
float xMin = 0.0;
float xMax = 800.0;
float yMin = 200.0;
float yMax = 600.0;

class Cell {
    PVector loc;
    float radius;
    int id;
    PVector v = new PVector(0,0,0);
    boolean migratory;

    Cell(float tempX, float tempY, float tempRad, int tempId,
	 boolean tempMigratory){

	loc = new PVector(tempX, tempY);
	radius = tempRad;
	id = tempId;
	migratory = tempMigratory;
    }
    void move(){
	v = loc.get();
	PVector repForce = sumFrep();
	repForce.mult(repConst);
	loc.add(repForce);
		PVector migForce = sumFmig();
	migForce.mult(migConst);
	loc.add(migForce);
	PVector fluForce = fFlu();
	fluForce.mult(fluConst);
	loc.add(fluForce);
	PVector adhForce = sumFadh();
	adhForce.mult(adhConst);
	loc.add(adhForce);
	float tx = constrain(loc.x,xMin,xMax);
	float ty = constrain(loc.y,yMin,yMax);
	loc.x = tx;
	loc.y = ty;
	v.sub(loc);
	v.mult(-1);
    }

    PVector sumFrep(){
	PVector total = new PVector(0,0,0);
	for(int i = 0; i < cells.size(); i++){
	    if(i != id){  //don't add self!
		Cell j = (Cell) cells.get(i);
		total.add(fRep(j));
	    }
	}
//	print(total);
	return total;
    }
    PVector fRep(Cell j){
	PVector jLoc = j.getLoc();
	PVector delta = new PVector (loc.x - j.getLoc().x, loc.y - j.getLoc().y);
	float overlap = (radius + j.getRadius()) - delta.mag();
	if (overlap > 0 && overlap != radius){//don't interact with self
	    delta.mult(k * (overlap)/(delta.mag() + epsilon));
	    return delta;
	}
	else
	    return new PVector(0.0,0.0,0.0);
    }

    PVector fAdh(Cell j){
	if(migratory && j.getMigratory()){
	    PVector jLoc = j.getLoc();
	    PVector delta = new PVector (loc.x - j.getLoc().x, 
					 loc.y - j.getLoc().y);
	    delta.mult(-1/(delta.mag() + epsilon));
	    return delta;
	}
	else
	    return new PVector(0.0,0.0,0.0);
    }

    PVector fMig(Cell j){
	if(migratory ^ j.getMigratory()){ //logical XOR!
	    PVector jLoc = j.getLoc();
	    PVector delta = new PVector (loc.x - j.getLoc().x, 
					 loc.y - j.getLoc().y);	    
	    float deltaNorm = delta.mag();
	    if(deltaNorm < 50){
		float dy = delta.y;
		float dx = delta.x;
		
		int s;
		if(migratory)
		    s = 1;
		else
		    s = - 1;
		return new PVector(s * abs(dy)/(deltaNorm + epsilon), 
				   s * (-abs(dy)/(dy + epsilon)) * 
				   (dx/(deltaNorm + epsilon)));
	    }
	    else
		return new PVector(0,0,0);
	}
	else
	    return new PVector(0,0,0);
    }
    
    PVector fFlu(){
	float u1 = random(0,1);
	float u2 = random(0,1);
	float k = sqrt(-2*log(u2));
	return new PVector(k * cos(u1), k * sin(u1));
    }

    PVector sumFmig(){
	PVector total = new PVector(0,0,0);
	for(int i = 0; i < cells.size(); i++){
	    if(i != id){  //don't add self!
		Cell j = (Cell) cells.get(i);
		total.add(fMig(j));
	    }
	}
//	print(total);
	return total;
    }

    PVector sumFadh(){
	PVector total = new PVector(0,0,0);
	for(int i = 0; i < cells.size(); i++){
	    if(i != id){  //don't add self!
		Cell j = (Cell) cells.get(i);
		total.add(fAdh(j));
	    }
	}
	return total;
    }

    void drawSelf(){
	if(migratory && id == 0)
	    fill(0,0,127);
	else if(migratory)
	    fill(127,0,0);
	else
	    fill(0,127,0);
	ellipse(loc.x,loc.y,radius * 2,radius * 2);
    }

    float getRadius(){return radius;}
    PVector getLoc()   {return loc;   }
    int getId()   {return id;   }
    boolean getMigratory()   {return migratory;   }
    PVector getV()   {return v;   }
}

void setup() 
{
  size(800, 800);
  //initialize cells array
  cells = new ArrayList();
  boolean mobile;
  PVector location,r;
  float disperse;
  for (int i = 0; i < numCells; i++){
      mobile = i % (numCells / numMigrants) == 0;
      location = mobile ? borderLocation : nurseLocation;	  
      disperse = mobile ? borderDisperse : nurseDisperse;	  
      r = new PVector(random(-disperse,disperse),
		      random(-disperse,disperse));
      r.add(location);
      cells.add(new Cell(r.x,r.y,radius,i,mobile));
  }

  for (int i = 0; i < cells.size(); i++){
      Cell c = (Cell) cells.get(i);
//      print("setting up cell #" + c.getId() + " i: " + i + "\n");
      c.drawSelf();

  }
}
int gen = 0;
Cell c;
PVector mu = new PVector(0,0,0);
PVector r;
PVector angMomentum = new PVector(0,0,0);
PVector linMomentum = new PVector(0,0,0);
float energy;
void draw() 
{
  background(51);
  //  print(gen + "\n");
  gen++;
  
  for (int i = 0; i < cells.size(); i++){
      c = (Cell) cells.get(i);
//      print("moving cell #" + c.getId() + " i: " + i + "\n");
      c.move();
      c.drawSelf();
//      print(c.getLoc());
      if(c.getMigratory()){
      mu.add(c.getLoc());
      }
  }
  mu.mult(1.0/numCells);
  angMomentum.set(0.0,0.0,0.0);
  linMomentum.set(0.0,0.0,0.0);
  energy = 0;
  for (int i = 0; i < cells.size(); i++){
      c = (Cell) cells.get(i);
      //      print("cell loc: " + c.getLoc()+"\n");
      //      print("cell velocity: " + c.getV()+"\n");
      if(c.getMigratory()){
      r = PVector.sub(c.getLoc(),mu);
      angMomentum.add(r.cross(c.getV()));
      linMomentum.add(c.getV());
      energy += pow((c.getV().mag()),2);
      }
  }
  //  print(momentum.z + "\n");
  println(linMomentum);
}

