nodes = []
distances = [] //2d array
newdistances = []
newnodes = []
newnodes2 = []
numberofnodes = 9
threshold = 10
offset = 1
canWidth = 700
canHeight = 700
speedupfactor = 25
repetitions = 1
speedofsound = 343
let pulserController
let stateController
let wind
let nodeCluster
let index1
noiseMag = 0
function setup() {
  createCanvas(canWidth, canHeight);
  seed = random(0, 1000)
  print("Seed: ", seed)
  randomSeed(seed)
  // randomSeed(180.0841701951379)
  // randomSeed(96.1659040443814)
  randomSeed(593.5639952957688)
  background(220);
  line(0, canHeight/2, canWidth, canHeight/2);
  line(canWidth/2, 0, canWidth/2, canHeight);
  stepOnePre();
  stepTwoPre();
  stepThreePre();
  stepFourPre();
  // tmpvec1 = new Vector(5.0, 3.0)
  // out = basisTransformation(new Vector(5.0, 3.0), tmpvec1, tmpvec1.eulerrotate(90))
  wind = new Vector(-1*noiseMag,noiseMag)
  // print("out: ", out[0],", ",out[1])
}


function draw() {
  index1 = int(max((mouseY-(canHeight/2))/((numberofnodes+1)*2.5), 0))
  index1 = int(min(index1, numberofnodes))
  // index = 7
  // print(index)
  testrows = iloc(newdistances[index1], 5)
  stepOne();
  stepTwo();
  stepThree();
  stepFour();
  push()
  stroke('green')
  fill('green')
  line(width/2, height/2, width/2 + 10*wind.x, height/2 + 10*wind.y)
  pop()
  // circle(mouseX, mouseY, 5)
}
/// given a matrix of distances between some points, returns the
/// point coordinates that best approximate the distances
function mds(distances, dimensions) {
    dimensions = dimensions || 2;

    // square distances
    var M = numeric.mul(-0.5, numeric.pow(distances, 2));

    // double centre the rows/columns
    function mean(A) { return numeric.div(numeric.add.apply(null, A), A.length); }
    var rowMeans = mean(M),
        colMeans = mean(numeric.transpose(M)),
        totalMean = mean(rowMeans);

    for (var i = 0; i < M.length; ++i) {
        for (var j =0; j < M[0].length; ++j) {
            M[i][j] += totalMean - rowMeans[i] - colMeans[j];
        }
    }

    // take the SVD of the double centred matrix, and return the
    // points from it
    var ret = numeric.svd(M),
        eigenValues = numeric.sqrt(ret.S);
    return ret.U.map(function(row) {
        return numeric.mul(row, eigenValues).splice(0, dimensions);
    })
}

function createOnes(dimensions){
  var ones = numeric.identity(dimensions)
    for(var i = 0; i < ones.length; ++i){
      for(var j = 0; j < ones[0].length; ++j){
        ones[i][j] = 1
      }
    }
  return ones
}

function cleanmds(distances, dimensions) {
    dimensions = dimensions || 2;

    // square distances
    var D2 = numeric.pow(distances, 2);

    // centering matrix C
    var ones = createOnes(numberofnodes+1)
    var C = numeric.sub(numeric.identity(numberofnodes+1),numeric.mul(1/(numberofnodes+1), ones))

    // Double centered B
    B = numeric.mul(-0.5, numeric.dot(C, numeric.dot(D2, C)))

    // take the SVD of the double centred matrix
    var ret = numeric.svd(B),
        eigenValues = numeric.sqrt(ret.S);
  
    // return points using SVD
    return ret.U.map(function(row) {
        return numeric.mul(row, eigenValues).splice(0, dimensions);
    })
}


class Vector{
 constructor(x, y){
   this.x = x
   this.y = y
 } 
  unitvec(){
    let magnitude = mag(this.x, this.y)
    return new Vector(this.x / magnitude, this.y / magnitude)
  }
  eulerrotate(nd){
    let rad = nd*PI/180;
    return new Vector(this.x*cos(rad)-this.y*sin(rad), this.y*cos(rad)+this.x*sin(rad))
  }
  add(vec2){
    return new Vector(this.x + vec2.x, this.y+vec2.y)
  }
  negative(){
    return new Vector(-1*this.x, -1*this.y)
  }
  magnitude(){
    return mag(this.x, this.y)
  }
  dot(vec2){
    return this.x*vec2.x+ this.y*vec2.y
  }
}

function basisTransformation(inputVec, basis1, basis2){
  //The idea is to transform the inputVec into a basis of {basis1, basis2},
  //here where basis1 is the vector of the wind and basis 2 is perpendicular to the wind
  // print(basis1)
  unitbasis1 = basis1.unitvec()
  unitbasis2 = basis2.unitvec()
  //Have to solve the system of linear equations ub1*i + ub2*j = iv
  //Or, ub1.x*i + ub2.x*j = iv.x, ub1.y*i + ub2.y*j = iv.y
  // j = (iv.x - ub1.x*i) /ub2.x
  // => ub1.y*i + ub2.y*(iv.x - ub1.x*i) / ub2.x = iv.y
  // => ub1.y*i + ub2.y*iv.x/ub2.x - ub2.y*ub1.x*i/ub2.x = iv.y
  // => (ub1.y-(ub2.y*ub1.x/ub2.x)) * i = iv.y - (ub2.y*iv.x/ub2.x)
  // => i = (iv.y - (ub2.y*iv.x/ub2.x)) / (ub1.y-(ub2.y*ub1.x/ub2.x))
  i = (inputVec.y - (unitbasis2.y * inputVec.x / unitbasis2.x)) / (unitbasis1.y - (unitbasis2.y*unitbasis1.x/unitbasis2.x))
  j = (inputVec.x - unitbasis1.x*i) / unitbasis2.x
  print("inputVec: ", inputVec.x, ", ", inputVec.y)
  print("Output: ", unitbasis1.x*i+unitbasis2.x*j, ", ",unitbasis2.y*i+unitbasis2.y*j)
  return [i, j]
  
}

function solveTwoCircleIntersection(x1, y1, x2, y2, r1, r2, R, sign){
  //   (x,y) = .5(x1+x2, y1+y2) + (r1^2 -r2^2)/(2R^2)(x2-x1, y2-y1) +- .5sqrt( 2*(r1^2+r2^2)/(R^2) - ((r1^2 - r2^2)^2)/R^4 -1)(y2-y1, x1-x2) where R = dist(c1, c2)
  x3 = 0.5*(x1+x2) + (r1*r1 - r2*r2)*(x2-x1)/(2*R*R) + sign*0.5*sqrt((2*(r1*r1+r2*r2)/(R*R)) - (pow(r1*r1 - r2*r2,2)/(R*R*R*R)) - 1)*(y2-y1)
  y3 = 0.5*(y1+y2) + (r1*r1 - r2*r2)*(y2-y1)/(2*R*R) + sign*0.5*sqrt((2*(r1*r1+r2*r2)/(R*R)) - (pow(r1*r1 - r2*r2,2)/(R*R*R*R)) - 1)*(x1-x2)
  //   print(0.5*(x1+x2) , (r1*r1 - r2*r2)*(x2-x1)/(2*R*R) ,(2*(r1*r1+r2*r2)/(R*R)) , (pow(r1*r1 - r2*r2,2)/(R*R*R*R)) ,-1,(y2-y1))
  // print(x3, y3)
  return [x3, y3]
}
class StateController{
  //The state controller needs to be able to keep track of what its currently doing,
  //and to move informaton between steps.
  
  //The steps of the state controller should be to send out a bluetooth pulse from base, wait until its finished, then send out an acoustic signal from node 0. Whenever this acoustic signal hits a node, that node should send a bluetooth message out to base. Once all signals are done, do the same with node 1
  constructor(){
    this.status = 0
    this.substep = 0
    this.curnode = 0
    this.pulsecontroller = new PulserController()
    this.currep = 0
  }
  check(){
    this.pulsecontroller.check()
    if(this.status == 0){
      this.subStepOne()
    }else if(this.status == 1){
      this.subStepTwo()
    }else if(this.status == 2){
      this.subStepThree()
    }
  }
  
  //2c / ((1/tij) + (1/tji))
  subStepOne(){
    if(this.pulsecontroller.pulses.length != 0){
      return //should wait for signals to be done
    }
    if(this.substep == 0){  
      this.pulsecontroller.addPulse(canWidth/4, canHeight/2 - 10, 'B', this.curnode)
      this.substep = 1
    }else if(this.substep == 1 && this.curnode < numberofnodes+1){
      this.pulsecontroller.addPulse(nodes[this.curnode].x, nodes[this.curnode].y, 'R', this.curnode)
      this.curnode += 1
      this.substep = 0
    }else if(this.currep == repetitions-1){ //once each node has sent out an acoustic pulse and you're done, move to next step!
      this.status = 1
    }else{
      this.curnode = 0
      this.substep = 0
      this.currep += 1
      speedupfactor += 5
      wind = new Vector(random(-1*noiseMag, noiseMag), random(-1*noiseMag, noiseMag))
    }
  }
//     subStepTwo(){
//         newnodes.push(new Node(canWidth/4, 3*canHeight/4, 0))
//   //say for the sake of making it look alike, we have a second node position (otherwise, rot indepdent?)
//   newnodes.push(new Node(nodes[1].x, nodes[1].y+canHeight/2, 1))
//   // newnodes.push(new Node(nodes[2].x, nodes[2].y+canHeight/2, 2))
//   // newnodes.push(new Node(canWidth/4 + distances[0][1], 3*canHeight/4, 1))
//   let answer
//   answer = solveTwoCircleIntersection(newnodes[0].x, newnodes[0].y, newnodes[1].x, newnodes[1].y, distances[0][2], distances[1][2], newnodes[0].distance(newnodes[1]), 1)
//   x3 = answer[0], y3=answer[1]
//   newnodes.push(new Node(x3, y3, 2))
//   for(let i = 3; i < nodes.length; i++){
//   answer = solveTwoCircleIntersection(newnodes[0].x, newnodes[0].y, newnodes[1].x, newnodes[1].y, distances[0][i], distances[1][i], newnodes[0].distance(newnodes[1]), 1)
//   x3 = answer[0], y3=answer[1]
//   if(abs(dist(newnodes[2].x, newnodes[2].y, x3, y3) - distances[2][i]) < threshold){
//     newnodes.push(new Node(x3, y3, i))
//   }else{
//   answer = solveTwoCircleIntersection(newnodes[0].x, newnodes[0].y, newnodes[1].x, newnodes[1].y, distances[0][i], distances[1][i], newnodes[0].distance(newnodes[1]), -1)
//     x3 = answer[0], y3=answer[1]
//     if(abs(dist(newnodes[2].x, newnodes[2].y, x3, y3) - distances[2][i]) < threshold){
//     newnodes.push(new Node(x3, y3, i))
//   }else{
//     print("Error!!!!")
//   }
//   }
//     this.status = 2
//   }
  
  
//   let tmpsum = 0
//   for(let i = 0; i < newnodes.length; i++){
//     for(let j = 0; j < newnodes.length; j++){
//       distance = newnodes[i].distance(newnodes[j])
//       tmpsum += abs(nodes[i].distance(nodes[j]) - distance) 
//     }
//   }
//   print("Error in recreation 1: ", tmpsum)
//       this.status = 2
//     }
    subStepTwo(){
      let newdistance
      for(let i = 0; i < distances.length; i++){
        for(let j = i; j < distances[0].length; j++){
          newdistance = 2*speedofsound/(((1/distances[i][j]) + (1/distances[j][i])))
          // newdistance = (distances[i][j] + distances[j][i]) / 2
          
          // newdistances[i][j] = (distances[i][j] + distances[j][i]) / 2
          // newdistances[j][i] = (distances[i][j] + distances[j][i]) / 2
          newdistances[i][j] = newdistance
          newdistances[j][i] = newdistance
          // print(newdistance, (distances[i][j] + distances[j][i]) / 2)
        }
      }
      this.status = 2
    }
    
//     subStepThree(){
//     newnodes2.push(new Node(3*canWidth/4, 3*canHeight/4, 0))
//   //say for the sake of making it look alike, we have a second node position (otherwise, rot indepdent?)
//   newnodes2.push(new Node(nodes[1].x+canWidth/2, nodes[1].y+canHeight/2, 1))
//   // newnodes2.push(new Node(nodes[2].x+canWidth/2, nodes[2].y+canHeight/2, 2))
//   // newnodes2.push(new Node(3*canWidth/4 + distances[0][1], 3*canHeight/4, 1))
  
//   let answer
  
//   // answer = solveTwoCircleIntersection(newnodes2[0].x, newnodes2[0].y, newnodes2[1].x, newnodes2[1].y, distances[0][2], distances[1][2], newnodes2[0].distance(newnodes2[1]), -1)
//   // x3 = answer[0], y3=answer[1]
//   // newnodes2.push(new Node(x3, y3, 2))
  
//   newnodes2.push(new Node(nodes[2].x+canWidth/2, nodes[2].y+canHeight/2, 2))
// for(let i = 3; i < nodes.length; i++){
    
  
//   answer = solveTwoCircleIntersection(newnodes2[0].x, newnodes2[0].y, newnodes2[1].x, newnodes2[1].y, newdistances[0][i], newdistances[1][i], newnodes2[0].distance(newnodes2[1]), -1)
//   x3 = answer[0], y3=answer[1]
//   if(abs(dist(newnodes2[2].x, newnodes2[2].y, x3, y3) - newdistances[2][i]) < threshold){
//     newnodes2.push(new Node(x3, y3, i))
//   }else{
//   answer = solveTwoCircleIntersection(newnodes2[0].x, newnodes2[0].y, newnodes2[1].x, newnodes2[1].y, newdistances[0][i], newdistances[1][i], newnodes2[0].distance(newnodes2[1]), 1)
//     x3 = answer[0], y3=answer[1]
//     if(abs(dist(newnodes2[2].x, newnodes2[2].y, x3, y3) - newdistances[2][i]) < threshold){
//     newnodes2.push(new Node(x3, y3, i))
//   }else{
//     print("Error!!!!")
//   }
//   }
//   }
//   let tmpsum = 0
//   for(let i = 0; i < newnodes2.length; i++){
//     for(let j = 0; j < newnodes2.length; j++){
//       distance = newnodes2[i].distance(newnodes2[j])      
//       tmpsum += abs(nodes[i].distance(nodes[j]) - distance) 
//       if(abs(nodes[i].distance(nodes[j]) - distance)  > 0.5){
//         print("i,", i, ", j,",j,": ",abs(nodes[i].distance(nodes[j]) - distance) )
//         line(newnodes2[i].x, newnodes2[i].y, newnodes2[j].x, newnodes2[j].y)
//       }
//     }
//   }
//   print("Error in recreation 2: ", tmpsum)
//   noLoop()
//   frameRate(0)
//   this.status = 3
//   }
    subStepThree(){
  let coords2 = cleanmds(newdistances, 2)
  // print("Ans 2: ", coords2)
  // let newnodes4 = []
  let correctionx = -1*coords2[0][0]
  let correctiony = -1*coords2[0][1]
  nodeCluster = new NodeCluster(3*canWidth/4, 3*canHeight/4)
  for(let i = 0; i < coords2.length; i++){
    // newnodes4.push(new Node(3*canWidth/4 + coords2[i][0] + correctionx, 3*canHeight/4 + coords2[i][1] + correctiony, i))
    nodeCluster.addNode(new Node(coords2[i][0] + correctionx, coords2[i][1] + correctiony, i))
  }
  // for(let i = 0; i < newnodes4.length; i++){
  //   fill('red')
  //   newnodes4[i].display()
  // }
  // fill('red')
  //RMSE
  let tmpsum = 0
  for(let i = 0; i < nodeCluster.arr.length; i++){
    for(let j = 0; j < nodeCluster.arr.length; j++){
      distance = nodeCluster.arr[i].distance(nodeCluster.arr[j])
      tmpsum += pow(nodes[i].distance(nodes[j]) - distance, 2) 
    }
  }
  tmpsum = tmpsum / (nodeCluster.arr.length * nodeCluster.arr.length)
  tmpsum = sqrt(tmpsum)
  //x2 = x1*cos(theta) - y1*sin(theta)
  //y2 = y1*cos(theta) + x1*sin(theta)
  // cos(theta) = (x2 + y1*sin(theta)) / x1
  // cos(theta) = (y2 - x1*sin(theta)) / y1
  // (x2*y1 + y1*y1*sin(theta)) = (y2*x1 - x1*x1*sin(theta))
  // (y1*y1 + x1*x1)*sin(theta) = (-1*x2*y1 + x1*y2)
  // theta = sin^-1((-1*x2*y1 + x1*y2) / (y1*y1 + x1*x1))
  
  // nodeCluster.mirror()
  let x1 = nodeCluster.arr[1].x - nodeCluster.arr[0].x
  let y1 = nodeCluster.arr[1].y - nodeCluster.arr[0].y
  let x2 = nodes[1].x - nodes[0].x
  let y2 = nodes[1].y - nodes[0].y
  // print(x1, y1, x2, y2, (-1*x2*y1 + x1*y2) / (y1*y1 + x1*x1))
  let theta = asin((-1*x2*y1 + x1*y2) / (y1*y1 + x1*x1))
  // print("theta: ", theta, " radians")
  // nodeCluster.rotateAll(theta)
  nodeCluster.rotateAll(45*PI/180)
  // nodeCluster.rotateAll(90*PI/180)
  nodeCluster.centerx += 15
  nodeCluster.display()
  
  print("Error in recreation 4, RMSE: ", tmpsum)
      this.status=3
    }
}
  
  
function heard(i, j, r, step){ //schedule a bluetooth pulse from j saying it heard i
  if(i != j && abs(nodes[i].distance(nodes[j]) - r/2) < step/2){
    stateController.pulsecontroller.addPulse(nodes[j].x, nodes[j].y, 'B', j)
    vec1 = new Vector(nodes[i].x - nodes[j].x, nodes[i].y - nodes[j].y)
    vec2 = new Vector(nodes[j].x - nodes[i].x, nodes[j].y - nodes[i].y)
    // distances[i][j] = nodes[i].distance(nodes[j])
    // distances[j][i] = nodes[j].distance(nodes[i])
    // distances[i][j] = vec1.add(wind).magnitude()
    // distances[j][i] = vec2.add(wind).magnitude()
    distances[i][j] += (vec1.magnitude() + vec1.unitvec().dot(wind) + random(-1*noiseMag, noiseMag)) / (speedofsound*repetitions*2)
    // distances[j][i] = (vec2.magnitude() + vec2.unitvec().dot(wind) + random(-1*noiseMag, noiseMag))/1
  }
}
class PulserController{
  constructor(){
    this.pulses = []
  }
  check(){
    for(let i = 0; i < this.pulses.length; i++){
      this.pulses[i].check()
    }
    this.pulses = this.pulses.filter(pulse => pulse.status == true)
  }
  addPulse(x, y, mode, node){
    this.pulses.push(new Pulser(x, y, node, heard))
    if(mode == 'R'){
      this.pulses[this.pulses.length-1].pulseR(x, y)
    }else{
      this.pulses[this.pulses.length-1].pulseB(x, y)
    }
  }
}
class Pulser{
    constructor(x, y, node, callback){
      this.maxradius = sqrt(2)*canWidth
      this.status = false
      this.x = x
      this.y = y
      this.r = 0
      this.step = 2*speedupfactor
      this.mode = 'R'
      this.node = node
      this.callback = callback
    }
    pulseR(x, y){ //acoustic signals are red
      this.mode = 'R'
      this.status = true
      this.x = x
      this.y = y
      this.r = this.r + this.step
      for(let i = 0; i < nodes.length; i++){
        this.callback(this.node, i, this.r, this.step)
      }
      push()
      stroke('red')
      noFill()
      circle(this.x, this.y, this.r)
      pop()
    }
    pulseB(x, y){ //bluetooth signals are blue
      this.mode = 'B'
      this.status = true
      this.x = x
      this.y = y
      this.r = this.r + this.step*3
      push()
      stroke('blue')
      noFill()
      circle(this.x, this.y, this.r)
      pop()
    }
    check(){
      if(this.r < this.maxradius && this.status){
        if(this.mode == 'R'){
          this.pulseR(this.x, this.y)
        }else{
          this.pulseB(this.x, this.y)
        }
      }else{
        this.r = 0
        this.status = false
      }
    }
}

class Node{
  constructor(x, y, id){
    this.x = x
    this.y = y
    this.id = id
  }
  display(){
    circle(this.x, this.y, 5)
    textSize(20)
    text(this.id, this.x, this.y)
  }
  distance(node2){
    return dist(this.x, this.y, node2.x, node2.y)
  }
  circ(r){
    circle(this.x, this.y, r)
  }
}

class NodeCluster{
  constructor(centerx, centery){
    this.centerx = centerx
    this.centery = centery
    this.arr = []
  }
  addNode(node){
    this.arr.push(node)
  }
  display(){
    for(let i = 0; i < this.arr.length; i++){
      circle(this.centerx + this.arr[i].x, this.centery + this.arr[i].y, 5)
      textSize(20)
      text(this.arr[i].id, this.centerx + this.arr[i].x, this.centery + this.arr[i].y)
    }
  }
  rotateAll(theta){
    for(let i = 0; i < this.arr.length; i++){
      let x = this.arr[i].x
      let y = this.arr[i].y
      this.arr[i].x = x*cos(theta) -y*sin(theta)
      this.arr[i].y = y*cos(theta) + x*sin(theta)
    }
  }
  mirror(){
    for(let i = 0; i < this.arr.length; i++){
      this.arr[i].x = -1*this.arr[i].x
    }
  }
}

function stepOnePre(){
  nodes.push(new Node(canWidth/4, canHeight/4, 0))
  // nodes.push(new Node(canWidth/4 + random(25, canHeight/4-25), canHeight/4, 1))
  wh = ceil(sqrt(numberofnodes))
  incr = (canWidth/2 - 50) / wh
  for(let i = 0; i < numberofnodes; i++){
    // nodes.push(new Node(random(25, canHeight/2 - 25), random(25, canHeight/2 - 25), i+1))
    nodes.push(new Node( 50 + incr*(i%wh) + random(-10, 10), 50 + incr*(floor((i/wh)%wh)) + random(-10,10) ,i+1))
  }
  fill(220)
  rect(0, 0, canWidth/2, canHeight/2)
  
  fill(0)
  nodes[0].display()
  fill(255)
  for(let i = 1; i < nodes.length; i++){
    nodes[i].display()
  }
  line(nodes[0].x, nodes[0].y, (nodes[numberofnodes].x - nodes[0].x)/4 + nodes[0].x, (nodes[numberofnodes].y - nodes[0].y)/4 + nodes[0].y)
  stateController = new StateController()
  stateController.check()
}
//First idea is in the top left quadrant, make n points randomy distributed around a central node
function stepOne(){
  
  fill(220)
  rect(0, 0, canWidth/2, canHeight/2)
  
//   fill(0)
//   nodes[0].display()
//   fill(255)
  for(let i = 0; i < nodes.length; i++){
    push()
    fill(0)
    if(testrows.includes(i)){
      // fill('red')
    }
    if(i == index1){
      // fill('blue')
    }
    nodes[i].display()
    pop()
  }
  // line(nodes[0].x, nodes[0].y, (nodes[numberofnodes].x - nodes[0].x)/4 + nodes[0].x, (nodes[numberofnodes].y - nodes[0].y)/4 + nodes[0].y)
  // node.display()
  stateController.check()
}
function stepTwoPre(){
  fill(220)
  rect(canWidth/2, 0, canWidth, canHeight/2)
  increment = ((canWidth/2)-5) / (nodes.length)
  fill(0)
  textSize(10)
  for(let i = 0; i < nodes.length; i++){
    tmp = []
    for(let j = 0; j < nodes.length; j++){
      // distance = nodes[i].distance(nodes[j])
      distance = 0
      tmp.push(0)
      // if(j >= i)
//         text(round(distance,3), canWidth/2 + 5 + increment*(j), 10 + increment*(i))
    }
    distances.push(tmp)
  }
}
function stepTwo(){
  fill(220)
  rect(canWidth/2, 0, canWidth, canHeight/2)
  increment = ((canWidth/2)-5) / (nodes.length)
  fill(0)
  textSize(10)
  for(let i = 0; i < nodes.length; i++){
    // tmp = []
    for(let j = 0; j < nodes.length; j++){
      // distance = nodes[i].distance(nodes[j])
      // tmp.push(distance)
      
      // if(j >= i)
        // text(round(distances[i][j],3), canWidth/2 + 5 + increment*(j), 10 + increment*(i))
    }
    // distances.push(tmp)
  }
}
// function stepThreePre(){
//   fill(220)
//   rect(0, canHeight/2, canWidth/2, canHeight)
// }
// function stepThree(){
//   fill(220)
//   rect(0, canHeight/2, canWidth/2, canHeight)
  
//   if(stateController.status > 1){
    
//   fill(0)
//   newnodes[0].display()
//   fill(255)
//   for(let i = 1; i < newnodes.length; i++){
//     newnodes[i].display()
//   }
//   }
// }
function stepThreePre(){
  fill(220)
  offsetx = 0
  offsety = canHeight/2
  offsetx2 = canHeight/2
  offsety2 = 0
  rect(offsetx, offsety, offsetx + canWidth/2, offsety + canHeight/2)
  increment = ((canWidth/2)-5) / (nodes.length)
  fill(0)
  textSize(10)
  for(let i = 0; i < nodes.length; i++){
    tmp = []
    for(let j = 0; j < nodes.length; j++){
      // distance = nodes[i].distance(nodes[j])
      distance = 0
      tmp.push(0)
      // if(j >= i)
        text(round(distance), offsetx2 + 5 + increment*(j), offsety2 + 10 + increment*(i))
    }
    newdistances.push(tmp)
  }
}
function indexofmin(array, cantbe){
  let smallest = 10000
  let smallestindex = -1
  for(let i = 0; i < array.length; i++){
    if(array[i] < smallest && !cantbe.includes(i)){
      smallestindex = i
      smallest = array[i]
    }
  }
  return smallestindex
}
function iloc(array, k){
  //array is a row of numbers. Return the indexes of the k smallest ones
  output = []
  modifiedarray = []
  for (i = 0; i < array.length; i++) {
  modifiedarray[i] = array[i];
}
  // print(modifiedarray)
  for(let i = 0; i < k; i++){
  index = indexofmin(modifiedarray, output)
  output.push(index)
  }
  // print(output)
  return output
}
function stepThree(){
  fill(220)
  offsetx = 0
  offsety = canHeight/2
  offsetx2 = canHeight/2
  offsety2 = 0
  rect(0, offsety, offsetx + canWidth/2, offsety + canHeight/2)
  increment = ((canWidth/2)-5) / (nodes.length)
  fill(0)
  textSize(10)
  
  
  for(let i = 0; i < nodes.length; i++){
    // tmp = []
    for(let j = 0; j < nodes.length; j++){
      // distance = nodes[i].distance(nodes[j])
      // tmp.push(distance)
      
      // if(j >= i)
        push()
        if(testrows.includes(i) && testrows.includes(j)){
          // fill('red')
        }
        if(i == index1 && j == index1){
          // print(index1)
          // fill('blue')
          // circle(offsetx + 5 + increment*j, offsety + 10 + increment*i, 5)
        }
        text(round(newdistances[i][j]), offsetx2 + 5 + increment*(j), offsety2 + 10 + increment*(i))
        pop()
    }
    // distances.push(tmp)
  }
}
function stepFourPre(){
  fill(220)
  rect(canWidth/2, canHeight/2, canWidth, canHeight)
  
  
}
function stepFour(){
  fill(220)
  rect(canWidth/2, canHeight/2, canWidth, canHeight)
  if(stateController.status > 2){
    fill(0)
  // newnodes2[0].display()
  // for(let i = 1; i < newnodes2.length; i++){
  //   newnodes2[i].display()
  // }
    // nodeCluster.display()
    for(let i = 0; i < nodes.length; i++){
    push()
    translate(nodeCluster.centerx, nodeCluster.centery)
    fill(0)
    if(testrows.includes(i)){
      // fill('red')
    }
    if(i == index1){
      // fill('blue')
    }
    // nodes[i].display()
    nodeCluster.arr[i].display()
    pop()
  }
  fill(255)
  }
}


class Matrice{
  constructor(h, w){
    this.mat = []
    this.h = h
    this.w = w
    for(let i = 0; i < this.h; i++){
      let tmpvec = []
      for(let j = 0; j < this.w; j++){
        tmpvec.push(0)
      }
      this.mat.push(tmpvec)
    }
  }
  populate(array){
    let count = 0
    if(array.length != this.h*this.w){
      print("Cannot populate matrice!")
      return -1
    }
    for(let i = 0; i < this.h; i++){
      for(let j = 0; j < this.w; j++){
        this.mat[i][j] = array[count]
        count++
      }
    }
  }
  populateMatrix(matrix){
    if(this.h != matrix.length || this.w != matrix[0].length){
      print("mismatch in matrice dimensions! Cannot populate")
    }
    for(let i = 0; i < this.h; i++){
      for(let j = 0; j < this.w; j++){
        this.mat[i][j] = matrix[i][j]
      }
    }
  }
  add(matrice2){
    if(this.h != matrice2.h || this.w != matrice2.w){
      print("Mismatch in matrice dimensions! Cannot add")
      return -1
    }
    for(let i = 0; i < this.h; i++){
      for(let j = 0; j < this.w; j++){
        this.mat[i][j] += matrice2[i][j]
      }
    }
  }
  subtract(matrice2){
    if(this.h != matrice2.h || this.w != matrice2.w){
      print("Mismatch in matrice dimensions! Cannot add")
      return -1
    }
    for(let i = 0; i < this.h; i++){
      for(let j = 0; j < this.w; j++){
        this.mat[i][j] -= matrice2[i][j]
      }
    }
  }
  addRow(row, index){
    if(this.w != row.length){
      print("Mismatch in row length and width!")
      return -1
    }
    for(let j = 0; j < this.w; j++){
      this.mat[index][j] += row[j]
    }
  }
  subtractRow(row, index){
    if(this.w != row.length){
      print("Mismatch in row length and width!")
      return -1
    }
    for(let j = 0; j < this.w; j++){
      this.mat[index][j] -= row[j]
    }
  }
  addColumn(column, index){
    if(this.h != column.length){
      print("Mismatch in column length and height!")
      return -1
    }
    for(let i = 0; i < this.h; i++){
      this.mat[i][index] += column[i]
    }
  }
  subtractColumn(column, index){
    if(this.h != column.length){
      print("Mismatch in column length and height!")
      return -1
    }
    for(let i = 0; i < this.h; i++){
      this.mat[i][index] -= column[i]
    }
  }
  multiply(constant){
    for(let i = 0; i < this.h; i++){
      for(let j = 0; j < this.w; j++){
        this.mat[i][j] *= constant
      }
    }
  }
  row(index){
    let tmpvec = []
    for(let j = 0; j < this.w; j++){
      tmpvec.push(this.mat[index][j])
    }
    return tmpvec
  }
  column(index){
    let tmpvec = []
    for(let i = 0; i < this.h; i++){
      tmpvec.push(this.mat[i][index])
    }
    return tmpvec
  }
  printout(){
    for(let i = 0; i < this.h; i++){
        print("[", this.mat[i], "]")
    }
  }
}