/* global noise */

// Set these two next lines to the ratio of your design and the prefix you want to use
const ratio = 1 / 1.4140// Canvas ratio (width / height, i.e. 3/4, 16/9, 1/1, 1/1.414 (Ax paper size))
const prefix = '81-82-83-84' // The filename we use when saving an image

const features = {} //  To hold the information we need to draw the design
let resizeTmr = null // This keeps track of the resize timer so we can cancel it if we need to
let thumbnailTaken = false // Keeps track of if we've taken a thumbnail yet

// Optional extras. The first two are grabbing the parameters from the URL
const urlSearchParams = new URLSearchParams(window.location.search)
const urlParams = Object.fromEntries(urlSearchParams.entries())
let forceDownloaded = false // Marked if we've downloaded the when told via the URL, so it only happens once

// This is also optional, and it used if we are animating the design
const animated = false
let nextFrame = null // requestAnimationFrame, and the ability to clear it

window.alba.params.seed = window.alba._testSeed()
// window.alba.params.seed = '0x181a992098e2ab4e00899b940cfb5b63dc67d154c5882b41b82bb038587f7760'
console.log(window.alba.params.seed)

class Rand {
  constructor (seed) {
    // PRNG from Piter
    const S = Uint32Array.of(9, 7, 5, 3)
    // eslint-disable-next-line no-return-assign
    this.prng = (a = 1) =>
      a *
      ((a = S[3]),
      (S[3] = S[2]),
      (S[2] = S[1]),
      (a ^= a << 11),
      (S[0] ^= a ^ (a >>> 8) ^ ((S[1] = S[0]) >>> 19)),
      S[0] / 2 ** 32);
    [...`${seed}`].map((c) =>
      this.prng((S[3] ^= c.charCodeAt() * 23205))
    )
  }

  r_d () {
    // random between 0 and 1
    return this.prng()
  }

  r_n (a, b) {
    // random float between a and b
    return a + (b - a) * this.r_d()
  }

  r_i (a, b) {
    // random int between a and b
    return ~~this.r_n(a, b + 1)
  }

  r_b (p) {
    // random boolean with probability of p
    return this.r_d() < p
  }

  r_c (list) {
    // random choice from list
    return list[this.r_i(0, list.length - 1)]
  }
}
const R = new Rand(window.alba.params.seed)

/*
 * A speed-improved perlin and simplex noise algorithms for 2D.
 *
 * Based on example code by Stefan Gustavson (stegu@itn.liu.se).
 * Optimisations by Peter Eastman (peastman@drizzle.stanford.edu).
 * Better rank ordering method by Stefan Gustavson in 2012.
 * Converted to Javascript by Joseph Gentle.
 *
 * Version 2012-03-09
 *
 * This code was placed in the public domain by its original author,
 * Stefan Gustavson. You may use it as you see fit, but
 * attribution is appreciated.
 *
 */
const module = window.noise = {}

function Grad (x, y, z) {
  this.x = x
  this.y = y
  this.z = z
}

Grad.prototype.dot2 = function (x, y) {
  return this.x * x + this.y * y
}

Grad.prototype.dot3 = function (x, y, z) {
  return this.x * x + this.y * y + this.z * z
}

const grad3 = [new Grad(1, 1, 0), new Grad(-1, 1, 0), new Grad(1, -1, 0), new Grad(-1, -1, 0),
  new Grad(1, 0, 1), new Grad(-1, 0, 1), new Grad(1, 0, -1), new Grad(-1, 0, -1),
  new Grad(0, 1, 1), new Grad(0, -1, 1), new Grad(0, 1, -1), new Grad(0, -1, -1)
]

const p = [151, 160, 137, 91, 90, 15,
  131, 13, 201, 95, 96, 53, 194, 233, 7, 225, 140, 36, 103, 30, 69, 142, 8, 99, 37, 240, 21,
  10,
  23,
  190, 6, 148, 247, 120, 234, 75, 0, 26, 197, 62, 94, 252, 219, 203, 117, 35, 11, 32, 57, 177,
  33,
  88, 237, 149, 56, 87, 174, 20, 125, 136, 171, 168, 68, 175, 74, 165, 71, 134, 139, 48, 27,
  166,
  77, 146, 158, 231, 83, 111, 229, 122, 60, 211, 133, 230, 220, 105, 92, 41, 55, 46, 245, 40,
  244,
  102, 143, 54, 65, 25, 63, 161, 1, 216, 80, 73, 209, 76, 132, 187, 208, 89, 18, 169, 200,
  196,
  135, 130, 116, 188, 159, 86, 164, 100, 109, 198, 173, 186, 3, 64, 52, 217, 226, 250, 124,
  123,
  5, 202, 38, 147, 118, 126, 255, 82, 85, 212, 207, 206, 59, 227, 47, 16, 58, 17, 182, 189,
  28,
  42,
  223, 183, 170, 213, 119, 248, 152, 2, 44, 154, 163, 70, 221, 153, 101, 155, 167, 43, 172, 9,
  129, 22, 39, 253, 19, 98, 108, 110, 79, 113, 224, 232, 178, 185, 112, 104, 218, 246, 97,
  228,
  251, 34, 242, 193, 238, 210, 144, 12, 191, 179, 162, 241, 81, 51, 145, 235, 249, 14, 239,
  107,
  49, 192, 214, 31, 181, 199, 106, 157, 184, 84, 204, 176, 115, 121, 50, 45, 127, 4, 150, 254,
  138, 236, 205, 93, 222, 114, 67, 29, 24, 72, 243, 141, 128, 195, 78, 66, 215, 61, 156, 180
]
// To remove the need for index wrapping, double the permutation table length
const perm = new Array(512)
const gradP = new Array(512)

// This isn't a very good seeding function, but it works ok. It supports 2^16
// different seed values. Write something better if you need more seeds.
module.seed = function (seed) {
  if (seed > 0 && seed < 1) {
    // Scale the seed out
    seed *= 65536
  }

  seed = Math.floor(seed)
  if (seed < 256) {
    seed |= seed << 8
  }

  for (let i = 0; i < 256; i++) {
    let v
    if (i & 1) {
      v = p[i] ^ (seed & 255)
    } else {
      v = p[i] ^ ((seed >> 8) & 255)
    }

    perm[i] = perm[i + 256] = v
    gradP[i] = gradP[i + 256] = grad3[v % 12]
  }
}

module.seed(0)

/*
  for(var i=0; i<256; i++) {
      perm[i] = perm[i + 256] = p[i];
      gradP[i] = gradP[i + 256] = grad3[perm[i] % 12];
  } */

// Skewing and unskewing factors for 2, 3, and 4 dimensions
const F2 = 0.5 * (Math.sqrt(3) - 1)
const G2 = (3 - Math.sqrt(3)) / 6

const F3 = 1 / 3
const G3 = 1 / 6

// 2D simplex noise
module.simplex2 = function (xin, yin) {
  let n0, n1, n2 // Noise contributions from the three corners
  // Skew the input space to determine which simplex cell we're in
  const s = (xin + yin) * F2 // Hairy factor for 2D
  let i = Math.floor(xin + s)
  let j = Math.floor(yin + s)
  const t = (i + j) * G2
  const x0 = xin - i + t // The x,y distances from the cell origin, unskewed.
  const y0 = yin - j + t
  // For the 2D case, the simplex shape is an equilateral triangle.
  // Determine which simplex we are in.
  let i1, j1 // Offsets for second (middle) corner of simplex in (i,j) coords
  if (x0 > y0) { // lower triangle, XY order: (0,0)->(1,0)->(1,1)
    i1 = 1
    j1 = 0
  } else { // upper triangle, YX order: (0,0)->(0,1)->(1,1)
    i1 = 0
    j1 = 1
  }
  // A step of (1,0) in (i,j) means a step of (1-c,-c) in (x,y), and
  // a step of (0,1) in (i,j) means a step of (-c,1-c) in (x,y), where
  // c = (3-sqrt(3))/6
  const x1 = x0 - i1 + G2 // Offsets for middle corner in (x,y) unskewed coords
  const y1 = y0 - j1 + G2
  const x2 = x0 - 1 + 2 * G2 // Offsets for last corner in (x,y) unskewed coords
  const y2 = y0 - 1 + 2 * G2
  // Work out the hashed gradient indices of the three simplex corners
  i &= 255
  j &= 255
  const gi0 = gradP[i + perm[j]]
  const gi1 = gradP[i + i1 + perm[j + j1]]
  const gi2 = gradP[i + 1 + perm[j + 1]]
  // Calculate the contribution from the three corners
  let t0 = 0.5 - x0 * x0 - y0 * y0
  if (t0 < 0) {
    n0 = 0
  } else {
    t0 *= t0
    n0 = t0 * t0 * gi0.dot2(x0, y0) // (x,y) of grad3 used for 2D gradient
  }
  let t1 = 0.5 - x1 * x1 - y1 * y1
  if (t1 < 0) {
    n1 = 0
  } else {
    t1 *= t1
    n1 = t1 * t1 * gi1.dot2(x1, y1)
  }
  let t2 = 0.5 - x2 * x2 - y2 * y2
  if (t2 < 0) {
    n2 = 0
  } else {
    t2 *= t2
    n2 = t2 * t2 * gi2.dot2(x2, y2)
  }
  // Add contributions from each corner to get the final noise value.
  // The result is scaled to return values in the interval [-1,1].
  return 70 * (n0 + n1 + n2)
}

// 3D simplex noise
module.simplex3 = function (xin, yin, zin) {
  let n0, n1, n2, n3 // Noise contributions from the four corners

  // Skew the input space to determine which simplex cell we're in
  const s = (xin + yin + zin) * F3 // Hairy factor for 2D
  let i = Math.floor(xin + s)
  let j = Math.floor(yin + s)
  let k = Math.floor(zin + s)

  const t = (i + j + k) * G3
  const x0 = xin - i + t // The x,y distances from the cell origin, unskewed.
  const y0 = yin - j + t
  const z0 = zin - k + t

  // For the 3D case, the simplex shape is a slightly irregular tetrahedron.
  // Determine which simplex we are in.
  let i1, j1, k1 // Offsets for second corner of simplex in (i,j,k) coords
  let i2, j2, k2 // Offsets for third corner of simplex in (i,j,k) coords
  if (x0 >= y0) {
    if (y0 >= z0) {
      i1 = 1
      j1 = 0
      k1 = 0
      i2 = 1
      j2 = 1
      k2 = 0
    } else if (x0 >= z0) {
      i1 = 1
      j1 = 0
      k1 = 0
      i2 = 1
      j2 = 0
      k2 = 1
    } else {
      i1 = 0
      j1 = 0
      k1 = 1
      i2 = 1
      j2 = 0
      k2 = 1
    }
  } else {
    if (y0 < z0) {
      i1 = 0
      j1 = 0
      k1 = 1
      i2 = 0
      j2 = 1
      k2 = 1
    } else if (x0 < z0) {
      i1 = 0
      j1 = 1
      k1 = 0
      i2 = 0
      j2 = 1
      k2 = 1
    } else {
      i1 = 0
      j1 = 1
      k1 = 0
      i2 = 1
      j2 = 1
      k2 = 0
    }
  }
  // A step of (1,0,0) in (i,j,k) means a step of (1-c,-c,-c) in (x,y,z),
  // a step of (0,1,0) in (i,j,k) means a step of (-c,1-c,-c) in (x,y,z), and
  // a step of (0,0,1) in (i,j,k) means a step of (-c,-c,1-c) in (x,y,z), where
  // c = 1/6.
  const x1 = x0 - i1 + G3 // Offsets for second corner
  const y1 = y0 - j1 + G3
  const z1 = z0 - k1 + G3

  const x2 = x0 - i2 + 2 * G3 // Offsets for third corner
  const y2 = y0 - j2 + 2 * G3
  const z2 = z0 - k2 + 2 * G3

  const x3 = x0 - 1 + 3 * G3 // Offsets for fourth corner
  const y3 = y0 - 1 + 3 * G3
  const z3 = z0 - 1 + 3 * G3

  // Work out the hashed gradient indices of the four simplex corners
  i &= 255
  j &= 255
  k &= 255
  const gi0 = gradP[i + perm[j + perm[k]]]
  const gi1 = gradP[i + i1 + perm[j + j1 + perm[k + k1]]]
  const gi2 = gradP[i + i2 + perm[j + j2 + perm[k + k2]]]
  const gi3 = gradP[i + 1 + perm[j + 1 + perm[k + 1]]]

  // Calculate the contribution from the four corners
  let t0 = 0.6 - x0 * x0 - y0 * y0 - z0 * z0
  if (t0 < 0) {
    n0 = 0
  } else {
    t0 *= t0
    n0 = t0 * t0 * gi0.dot3(x0, y0, z0) // (x,y) of grad3 used for 2D gradient
  }
  let t1 = 0.6 - x1 * x1 - y1 * y1 - z1 * z1
  if (t1 < 0) {
    n1 = 0
  } else {
    t1 *= t1
    n1 = t1 * t1 * gi1.dot3(x1, y1, z1)
  }
  let t2 = 0.6 - x2 * x2 - y2 * y2 - z2 * z2
  if (t2 < 0) {
    n2 = 0
  } else {
    t2 *= t2
    n2 = t2 * t2 * gi2.dot3(x2, y2, z2)
  }
  let t3 = 0.6 - x3 * x3 - y3 * y3 - z3 * z3
  if (t3 < 0) {
    n3 = 0
  } else {
    t3 *= t3
    n3 = t3 * t3 * gi3.dot3(x3, y3, z3)
  }
  // Add contributions from each corner to get the final noise value.
  // The result is scaled to return values in the interval [-1,1].
  return 32 * (n0 + n1 + n2 + n3)
}

// ##### Perlin noise stuff

function fade (t) {
  return t * t * t * (t * (t * 6 - 15) + 10)
}

function lerp (a, b, t) {
  return (1 - t) * a + t * b
}

// 2D Perlin Noise
module.perlin2 = function (x, y) {
  // Find unit grid cell containing point
  let X = Math.floor(x)
  let Y = Math.floor(y)
  // Get relative xy coordinates of point within that cell
  x = x - X
  y = y - Y
  // Wrap the integer cells at 255 (smaller integer period can be introduced here)
  X = X & 255
  Y = Y & 255

  // Calculate noise contributions from each of the four corners
  const n00 = gradP[X + perm[Y]].dot2(x, y)
  const n01 = gradP[X + perm[Y + 1]].dot2(x, y - 1)
  const n10 = gradP[X + 1 + perm[Y]].dot2(x - 1, y)
  const n11 = gradP[X + 1 + perm[Y + 1]].dot2(x - 1, y - 1)

  // Compute the fade curve value for x
  const u = fade(x)

  // Interpolate the four results
  return lerp(
    lerp(n00, n10, u),
    lerp(n01, n11, u),
    fade(y))
}

// 3D Perlin Noise
module.perlin3 = function (x, y, z) {
  // Find unit grid cell containing point
  let X = Math.floor(x)
  let Y = Math.floor(y)
  let Z = Math.floor(z)
  // Get relative xyz coordinates of point within that cell
  x = x - X
  y = y - Y
  z = z - Z
  // Wrap the integer cells at 255 (smaller integer period can be introduced here)
  X = X & 255
  Y = Y & 255
  Z = Z & 255

  // Calculate noise contributions from each of the eight corners
  const n000 = gradP[X + perm[Y + perm[Z]]].dot3(x, y, z)
  const n001 = gradP[X + perm[Y + perm[Z + 1]]].dot3(x, y, z - 1)
  const n010 = gradP[X + perm[Y + 1 + perm[Z]]].dot3(x, y - 1, z)
  const n011 = gradP[X + perm[Y + 1 + perm[Z + 1]]].dot3(x, y - 1, z - 1)
  const n100 = gradP[X + 1 + perm[Y + perm[Z]]].dot3(x - 1, y, z)
  const n101 = gradP[X + 1 + perm[Y + perm[Z + 1]]].dot3(x - 1, y, z - 1)
  const n110 = gradP[X + 1 + perm[Y + 1 + perm[Z]]].dot3(x - 1, y - 1, z)
  const n111 = gradP[X + 1 + perm[Y + 1 + perm[Z + 1]]].dot3(x - 1, y - 1, z - 1)

  // Compute the fade curve value for x, y, z
  const u = fade(x)
  const v = fade(y)
  const w = fade(z)

  // Interpolate
  return lerp(
    lerp(
      lerp(n000, n100, u),
      lerp(n001, n101, u), w),
    lerp(
      lerp(n010, n110, u),
      lerp(n011, n111, u), w),
    v)
}

// A function that finds the intersection point of two lines
const faultyLineIntersectsEdge = (line1, line2) => {
  // Get the points of the lines
  const p1 = line1.p1
  const p2 = line1.p2
  const p3 = line2.p1
  const p4 = line2.p2
  // Get the x and y values of the lines
  const x1 = p1.x
  const x2 = p2.x
  const x3 = p3.x
  const x4 = p4.x
  const y1 = p1.y
  const y2 = p2.y
  const y3 = p3.y
  const y4 = p4.y
  // Calculate the denominator
  const denominator = ((x1 - x2) * (y3 - y4)) - ((y1 - y2) * (x3 - x4))
  // If the denominator is 0, the lines are parallel
  if (denominator === 0) {
    return null
  }
  // Calculate the intersection point
  const x = (((x1 * y2) - (y1 * x2)) * (x3 - x4) - ((x1 - x2) * (x3 * y4 - y3 * x4))) / denominator
  const y = (((x1 * y2) - (y1 * x2)) * (y3 - y4) - ((y1 - y2) * (x3 * y4 - y3 * x4))) / denominator
  // Check if the intersection point is on the line segments
  if (x < Math.min(x1, x2) || x > Math.max(x1, x2) || x < Math.min(x3, x4) || x > Math.max(x3, x4)) {
    return null
  }
  if (y < Math.min(y1, y2) || y > Math.max(y1, y2) || y < Math.min(y3, y4) || y > Math.max(y3, y4)) {
    return null
  }
  // Return the point of intersection
  return {
    x,
    y
  }
}

// We need a function that calculates if a line intersects a polygon
const faultyLineIntersectsPoly = (line, poly) => {
  let intersection1 = null
  let intersection2 = null
  // Loop through the points of the poly grabbing the edges
  for (let i = 0; i < poly.length; i++) {
    const p1 = poly[i]
    const p2 = (i < poly.length - 1) ? poly[i + 1] : poly[0]
    const edge = {
      p1: {
        x: p1[0],
        y: p1[1]
      },
      p2: {
        x: p2[0],
        y: p2[1]
      }
    }
    // Get the intersection point of the line and the edge
    const intersection = faultyLineIntersectsEdge(line, edge)
    // If there is an intersection, set it
    if (intersection) {
      if (!intersection1) {
        intersection1 = intersection
      } else {
        intersection2 = intersection
      }
    }
  }
  if (intersection1 && intersection2) {
    return {
      p1: intersection1,
      p2: intersection2
    }
  }
  return null
}

// This decides if we're going to keep the square or subdivide it
const subDivideSquare = (corners, depth, layer) => {
  // Work out the middle of the square
  const middle = {
    x: corners.tl.x + ((corners.tr.x - corners.tl.x) / 2),
    y: corners.tl.y + ((corners.bl.y - corners.tl.y) / 2)
  }
  // Grab a value from the perlin noise for the middle of the square
  const perlinValue = (1 + noise.perlin2(features.perlinOffsets[layer].x + (middle.x * features.perlinOffsets[layer].scale), features.perlinOffsets[layer].y + (middle.y * features.perlinOffsets[layer].scale))) / 2

  // Work out the chance of it being subdivided
  let subChance = 0.9
  for (let i = 0; i < depth; i++) subChance *= 0.666
  // Don't let his carry on forever
  if (subChance < 0.2) subChance = -1

  // Now work out if we're going to subdivide it
  if (perlinValue < subChance) {
    // We are going to subdivide it, first work out the length of half the square
    const halfLength = (corners.tr.x - corners.tl.x) / 2
    // Now do two loops to work out the four corners of the new squares
    for (let y = 0; y < 2; y++) {
      for (let x = 0; x < 2; x++) {
        const newCorners = {
          tl: {
            x: corners.tl.x + (x * halfLength),
            y: corners.tl.y + (y * halfLength)
          },
          tr: {
            x: corners.tl.x + ((x + 1) * halfLength),
            y: corners.tl.y + (y * halfLength)
          },
          bl: {
            x: corners.tl.x + (x * halfLength),
            y: corners.tl.y + ((y + 1) * halfLength)
          },
          br: {
            x: corners.tl.x + ((x + 1) * halfLength),
            y: corners.tl.y + ((y + 1) * halfLength)
          }
        }
        subDivideSquare(newCorners, depth + 1, layer)
      }
    }
  } else {
    const newSquare = {
      corners
    }
    newSquare.middle = {
      x: corners.tl.x + ((corners.tr.x - corners.tl.x) / 2),
      y: corners.tl.y + ((corners.bl.y - corners.tl.y) / 2)
    }
    newSquare.rotate = 180 * noise.perlin2(features.perlinRotationOffsets[layer].x + (newSquare.middle.x * features.perlinRotationOffsets[layer].scale), features.perlinRotationOffsets[layer].y + (newSquare.middle.y * features.perlinRotationOffsets[layer].scale))
    features.squares[layer].push(newSquare)
  }
}

// Now we want a function that picks an item from a weighted list
// the list is an array of items, where each item is an array of two values
// the first is what we want back the second is the weighting
const pickWeighted = (list) => {
  // First work out the total weighting
  let totalWeighting = 0
  for (let i = 0; i < list.length; i++) {
    totalWeighting += list[i][1]
  }
  // Now pick a random number between 0 and the total
  let randomNumber = R.prng() * totalWeighting
  // Now loop through the list again, subtracting the weighting each time
  // until we get to 0 or less
  let itemIndex = 0
  while (randomNumber > 0) {
    randomNumber -= list[itemIndex][1]
    itemIndex++
  }
  // Now we have the item index, we can grab the item
  return list[itemIndex - 1][0]
}

const palettes = [
  [{
    name: 'Ode to White Mono',
    primaryColours: [['#f2f2f2', 1]],
    backgroundColours: [['#2E5B7A', 0.1], ['#24A385', 0.1], ['#C22127', 0.1], ['#DB7630', 0.1], ['#1F1610', 0.1], ['#1A2F4B', 0.1], ['#F9D166', 0.1]]
  }, 0.09],

  [{
    name: 'Oil',
    primaryColours: [['#75FAC8', 0.2], ['#567C6A', 0.2], ['#92DEFF', 0.2], ['#7A96D4', 0.2], ['#E49BFD', 0.1], ['#655A63', 0.1]],
    backgroundColours: [['#000000', 0.1]]
  }, 0.043],
  [{
    name: 'Taurine',
    primaryColours: [['#A60000', 0.1], ['#CF0202', 0.1], ['#BE3C35', 0.1], ['#4A0A0B', 0.1], ['#DA9F2B', 0.2]],
    backgroundColours: [['#2B0B0B', 0.0666]]
  }, 0.03],
  [{
    name: 'Superdeluxe',
    primaryColours: [['#24A385', 0.1], ['#F9D166', 0.1], ['#DB7630', 0.1], ['#B1CFC9', 0.1], ['#F2A497', 0.1]],
    backgroundColours: [['#111B2F', 0.1]]
  }, 0.043],
  [{
    name: 'Shibuya',
    primaryColours: [['#CEA897', 0.1], ['#ffd700', 0.1], ['#E52225', 0.1], ['#e7007c', 0.1], ['#0057A9', 0.1], ['#5DD7F0', 0.1], ['#f9f2e0', 0.1]],
    backgroundColours: [['#1D1D1D', 0.1]]
  }, 0.043],
  [{
    name: 'Melody',
    primaryColours: [['#D42C2A', 0.1], ['#f9f2e0', 0.1], ['#F0AE02', 0.1], ['#0a7a41', 0.1], ['#2A529E', 0.1], ['#f7c1b4', 0.1], ['#5DD7F0', 0.1]],
    backgroundColours: [['#1d1d1c', 0.3], ['#F9F2E0', 0.1]]
  }, 0.043],
  [{
    name: 'Old Town Gas Station',
    primaryColours: [['#06C6EE', 0.3], ['#f0d31c', 0.3], ['#BF9922', 0.1], ['#f1603a', 0.3], ['#f1f1f2', 0.1]],
    backgroundColours: [['#3D2722', 0.1]]
  }, 0.043],
  [{
    name: 'Matrix',
    primaryColours: [['#1F4D2B', 0.1], ['#24BF57', 0.1], ['#81D090', 0.1], ['#93E6A3', 0.1]],
    backgroundColours: [['#0E0D1D', 0.1]]
  }, 0.03],

  [{
    name: 'Shibuya Reverse',
    primaryColours: [['#CEA897', 0.1], ['#ffd700', 0.1], ['#E52225', 0.1], ['#e7007c', 0.1], ['#0057A9', 0.1], ['#5DD7F0', 0.1], ['#1D1D1D', 0.1]],
    backgroundColours: [['#f9f2e0', 0.1]]
  }, 0.0666],
  [{
    name: 'Print Me',
    primaryColours: [['#ffd700', 0.1], ['#e7007c', 0.1], ['#5DD7F0', 0.1], ['#1D1D1D', 0.1]],
    backgroundColours: [['#f2f2f2', 0.1]]
  }, 0.0666],
  [{
    name: 'Fundamental Riso Shift',
    primaryColours: [['#2659A3', 0.1], ['#f6c0d4', 0.1], ['#E73B1D', 0.1], ['#ffcf10', 0.1], ['#1D1D1D', 0.1]],
    backgroundColours: [['#f9f2e0', 0.1], ['#F8F9F4', 0.1]]
  }, 0.0666],
  [{
    name: 'Notaiesque',
    primaryColours: [['#074690', 0.1], ['#1A80C0', 0.1], ['#70aee0', 0.1], ['#9CD2DC', 0.1], ['#f29ca7', 0.1]],
    backgroundColours: [['#F8F9F4', 0.1]]
  }, 0.0666],
  [{
    name: 'Earthen',
    primaryColours: [['#918679', 0.1], ['#d1c9be', 0.1], ['#594B42', 0.1], ['#d3a679', 0.1]],
    backgroundColours: [['#FAFAFA', 0.1]]
  }, 0.0666],
  [{
    name: 'Unshifted',
    primaryColours: [['#F15060', 0.2], ['#FFe800', 0.1], ['#00a95c', 0.1], ['#62c2b1', 0.1], ['#3d5588', 0.1], ['#ffb511', 0.1], ['#0078bf', 0.1]],
    backgroundColours: [['#E7CACF', 0.1]]
  }, 0.0666],
  [{
    name: 'Hanamura',
    primaryColours: [['#FF98B4', 0.3], ['#FF5735', 0.3], ['#FFCE49', 0.3], ['#EDE8DC', 0.1]],
    backgroundColours: [['#F7F4ED', 0.3], ['#0E0D1B', 0.1]]
  }, 0.0666],
  [{
    name: 'Drift',
    primaryColours: [['#4D6BAC', 0.1], ['#9B69A5', 0.1], ['#E52F36', 0.1], ['#F68C1A', 0.1]],
    backgroundColours: [['#F9F0DE', 0.1]]
  }, 0.0666],
  [{
    name: 'Breaker',
    primaryColours: [['#4D6BAC', 0.1], ['#487A9C', 0.1], ['#FDFDFD', 0.1], ['#EAF5E8', 0.1], ['#031A29', 0.1], ['#253C54', 0.1], ['#E51531', 0.3]],
    backgroundColours: [['#F9F0DE', 0.3], ['#000030', 0.1]]
  }, 0.0666]
]

/*
 * This is your setup function, it gets called right at the start, and only once.
 * This is where you make all your decisions about what the design is going to look like.
 * All the random choices should be made here, and the information stored in the features object
 * so that the drawCanvas() function can use it to draw the design.
 *
 * As you want to do more complicated things you'll want to move beyond this simple setup function,
 * but for the moment this is all we need (we'll cover more in a future YouTube video)
 *
 * The features object is global, so you can access it from anywhere in your code.
 */
const setup = () => {
  // We want to pick one of the palettes at random, but based on the weightings
  // So first we need to work out the total weightings
  let totalWeighting = 0
  for (let i = 0; i < palettes.length; i++) {
    totalWeighting += palettes[i].weighting
  }
  // Now pick a random number between 0 and the total
  let randomNumber = R.prng() * totalWeighting
  // Now loop through the palettes again, subtracting the weighting each time
  // until we get to 0 or less
  let paletteIndex = 0
  while (randomNumber > 0) {
    randomNumber -= palettes[paletteIndex].weighting
    paletteIndex++
  }
  // Now we have the palette index, we can grab the palette
  features.palette = pickWeighted(palettes)

  // Now we want to pick a background colour
  features.backgroundColour = pickWeighted(features.palette.backgroundColours)
  // Work out if on average the background is dark or light
  // split the background into RGB values from the hex
  const r = parseInt(features.backgroundColour.substring(1, 3), 16)
  const g = parseInt(features.backgroundColour.substring(3, 5), 16)
  const b = parseInt(features.backgroundColour.substring(5, 7), 16)
  // Work out the average
  const average = (r + g + b) / 3
  // If the average is less than 64
  features.textureColour = 'black'
  if (average < 128) features.textureColour = 'white'

  // Remove the background colour from the the primary colours if it's in there
  const newPrimaryColours = []
  for (let i = 0; i < features.palette.primaryColours.length; i++) {
    if (features.palette.primaryColours[i][0] !== features.backgroundColour) {
      newPrimaryColours.push(features.palette.primaryColours[i])
    }
  }
  features.palette.primaryColours = newPrimaryColours

  // Work out how croppy we are going to be with the lines
  if (R.prng() < 0.25) {
    features.cropLines = 0.8
    if (R.prng() < 0.0) features.cropLines = 1.0
  }
  // PART ONE, decide _ALL THE THINGS_ and put what we'll need access to
  // in the features object
  let percentOfLines = 1
  while (percentOfLines > 0.21) {
    features.layers = Math.floor((R.prng() * 6 + R.prng() * 6) / 2) + 2
    // If we have the print-me palette, we need at least 4 layers
    if (features.palette.name === 'Print Me' && features.layers < 4) features.layers = 4
    if (features.palette.name === 'Earthenware' && features.layers < 4) features.layers = 4

    // features.layers = 1
    // We need something to hold all the lines
    features.lineHolder = []
    features.linesSet = false

    // We are going to have some perlin noise so we'll have some offsets and scale
    features.perlinOffsets = []
    features.perlinRotationOffsets = []
    features.perlinCrop = []

    const pickedColours = []

    for (let i = 0; i < features.layers; i++) {
      // If the pickedColours length is less than the primary colours, then we
      // need to keep picking colours until we have one that isn't already in the list
      // otherwise we just pick one at random
      let colour = pickWeighted(features.palette.primaryColours)
      if (pickedColours.length < features.palette.primaryColours.length) {
        while (pickedColours.indexOf(colour) !== -1) {
          colour = pickWeighted(features.palette.primaryColours)
        }
      }
      pickedColours.push(colour)

      //   Push the lines into the line holder
      features.lineHolder.push({
        lines: [],
        colour
      })

      // If we are on the first layer, make it random, otherwise copy the values over
      if (i === 0) {
        features.perlinOffsets.push({
          x: R.prng() * 1000 + 4000,
          y: R.prng() * 1000 + 4000,
          scale: 8
        })

        features.perlinRotationOffsets.push({
          x: R.prng() * 1000 + 4000,
          y: R.prng() * 1000 + 4000,
          scale: 2
        })
        features.perlinCrop.push({
          x: R.prng() * 1000 + 4000,
          y: R.prng() * 1000 + 4000,
          scale: 8
        })
      } else {
      // Copy them over
        features.perlinOffsets.push({
          x: features.perlinOffsets[0].x + (i * 0.25),
          y: features.perlinOffsets[0].y + (i * 0.125),
          scale: features.perlinOffsets[0].scale + (i * 0.0)
        })
        features.perlinRotationOffsets.push({
          x: features.perlinRotationOffsets[0].x + (i * 0.01),
          y: features.perlinRotationOffsets[0].y + (i * 0.005),
          scale: features.perlinRotationOffsets[0].scale + (i * 0.01)
        })
        features.perlinCrop.push({
          x: features.perlinCrop[0].x + (i * 0.01),
          y: features.perlinCrop[0].y + (i * 0.005),
          scale: features.perlinCrop[0].scale + (i * 0.01)
        })
      }
    }

    // We are going to be using a 1 by 1 * ratio co-ordinate system, so we need to work within that
    // I want to work out a grid that nicely fills the page, making sure to have a border
    features.sideBorder = 0.05
    // Pick a number of squares across the page, somewhere between 4 and 7
    features.squaresAcross = pickWeighted([[3, 1], [5, 2], [6, 3], [7, 4], [10, 5]])
    // Work out the size of the squares
    features.squareSize = (1 - (features.sideBorder * 2)) / features.squaresAcross
    // Work out the number of squares down the page
    features.squaresDown = Math.floor(((1 / ratio) - (features.sideBorder * 2)) / features.squareSize)
    // Work out the topBottomBorder
    features.topBottomBorder = features.sideBorder

    // Swap squares across and down, now
    const temp = features.squaresAcross
    features.squaresAcross = features.squaresDown
    features.squaresDown = temp - Math.floor(features.squaresAcross / 8)

    // Now have an array of all the squares
    const calculateCorner = (x, y, size, borderX, borderY) => {
      return {
        x: x * size + borderX,
        y: y * size + borderY
      }
    }
    features.squares = []
    for (let i = 0; i < features.layers; i++) {
      features.squares.push([])
      for (let y = 0; y < features.squaresDown; y++) {
        for (let x = 0; x < features.squaresAcross; x++) {
        // Work out the four corners of the square
          const squareSize = features.squareSize
          const sideBorder = features.sideBorder
          const topBottomBorder = features.topBottomBorder

          const tl = calculateCorner(x, y, squareSize, sideBorder, topBottomBorder)
          const tr = calculateCorner(x + 1, y, squareSize, sideBorder, topBottomBorder)
          const bl = calculateCorner(x, y + 1, squareSize, sideBorder, topBottomBorder)
          const br = calculateCorner(x + 1, y + 1, squareSize, sideBorder, topBottomBorder)

          const corners = { tl, tr, bl, br }
          subDivideSquare(corners, 0, i)
        }
      }
    }

    // Now go through all the squares logging them to the console
    const outerSizeMod = 3
    const outerLines = 30
    let maxTotalSquares = 0
    let totalFinalLines = 0

    for (let layer = 0; layer < features.layers; layer++) {
      maxTotalSquares += features.squares[layer].length
      for (let i = 0; i < features.squares[layer].length; i++) {
        const square = features.squares[layer][i]
        const squareSize = square.corners.tr.x - square.corners.tl.x
        const squareRotCos = Math.cos(square.rotate * (Math.PI / 180))
        const squareRotSin = Math.sin(square.rotate * (Math.PI / 180))

        // Make a culling poly from the square
        const cullingPoly = [
          [square.corners.tl.x - square.middle.x, square.corners.tl.y - square.middle.y],
          [square.corners.tr.x - square.middle.x, square.corners.tr.y - square.middle.y],
          [square.corners.br.x - square.middle.x, square.corners.br.y - square.middle.y],
          [square.corners.bl.x - square.middle.x, square.corners.bl.y - square.middle.y]
        ]

        // Now we want to make a bunch of lines that fill a square bigger than the square
        const biggerSquareSize = (cullingPoly[1][0] - cullingPoly[0][0]) * outerSizeMod
        const halfBiggerSquareSize = biggerSquareSize / 2
        const lineStep = biggerSquareSize / outerLines
        // Make a bunch of lines
        const lines = []
        for (let i = 0; i < outerLines; i++) {
          const baseX = (-halfBiggerSquareSize) + (i * lineStep)

          const p1 = {
            x: baseX * squareRotCos + halfBiggerSquareSize * squareRotSin,
            y: -halfBiggerSquareSize * squareRotCos + baseX * squareRotSin
          }
          const p2 = {
            x: baseX * squareRotCos - halfBiggerSquareSize * squareRotSin,
            y: halfBiggerSquareSize * squareRotCos + baseX * squareRotSin
          }

          lines.push({ p1, p2 })
        }

        // Now move all the points in the lines to the middle of the square
        const maxBorderX = 1 - features.sideBorder
        const halfSquareSize = squareSize * 0.5
        let flipXY = false
        const newLines = lines.map(line => {
          const newLine = faultyLineIntersectsPoly(line, cullingPoly)
          if (!newLine) return null

          // Adjust line points and calculate midpoint
          const adjustPoint = (point) => ({
            x: point.x + square.middle.x,
            y: point.y + square.middle.y
          })
          newLine.p1 = adjustPoint(newLine.p1)
          newLine.p2 = adjustPoint(newLine.p2)
          const midPoint = {
            x: (newLine.p1.x + newLine.p2.x) / 2,
            y: (newLine.p1.y + newLine.p2.y) / 2
          }

          // Calculate new end points
          // But only extend them if the random number is greater than the cropLines value
          // Work out the chance of it being cropped based on the perlinCrop value
          let cropChance = (1 + noise.perlin2(features.perlinCrop[layer].x + (midPoint.x * features.perlinCrop[layer].scale), features.perlinCrop[layer].y + (midPoint.y * features.perlinCrop[layer].scale))) / 2
          if (features.cropLines) cropChance = features.cropLines
          if (R.prng() > cropChance) {
            const gradient = (newLine.p2.y - newLine.p1.y) / (newLine.p2.x - newLine.p1.x)
            const angle = Math.atan(gradient)
            const deltaX = Math.cos(angle) * halfSquareSize
            const deltaY = Math.sin(angle) * halfSquareSize
            newLine.p1 = { x: midPoint.x - deltaX, y: midPoint.y - deltaY }
            newLine.p2 = { x: midPoint.x + deltaX, y: midPoint.y + deltaY }
          }
          // Flip X and Y if necessary
          if (newLine.p1.x > maxBorderX || newLine.p2.x > maxBorderX) {
            flipXY = true;
            [newLine.p1.x, newLine.p1.y, newLine.p2.x, newLine.p2.y] =
              [newLine.p1.y, newLine.p1.x, newLine.p2.y, newLine.p2.x]
          }

          // Recalculate the middle of the line and add it to the line object
          newLine.middle = {
            x: (newLine.p1.x + newLine.p2.x) / 2,
            y: (newLine.p1.y + newLine.p2.y) / 2
          }

          return newLine
        }).filter(line => line !== null)

        // Add the lines to the square
        totalFinalLines += newLines.length
        square.lines = newLines
        // If we flipped the x and y values of the lines, we need to flip the x and y values of the
        // corners and middle of the square too
        if (flipXY) {
          [square.corners.tl.x, square.corners.tl.y] = [square.corners.tl.y, square.corners.tl.x];
          [square.corners.tr.x, square.corners.tr.y] = [square.corners.tr.y, square.corners.tr.x];
          [square.corners.bl.x, square.corners.bl.y] = [square.corners.bl.y, square.corners.bl.x];
          [square.corners.br.x, square.corners.br.y] = [square.corners.br.y, square.corners.br.x];
          [square.middle.x, square.middle.y] = [square.middle.y, square.middle.x]
        }
      }
    }
    // In theory the maximun number of lines is the number of squares * the number of lines
    // per square, but we're going to have less than that because of culling. Work out the percent
    // we have ended up with
    percentOfLines = totalFinalLines / (maxTotalSquares * outerLines)
  }

  // Now that we have all the squares and lines, I want to work out where the lowest y position
  // (which is really the highest value because of the way the co-ordinate system works)
  // of the middle of the squares is, so we can report it on the console.
  let lowestY = 0
  for (let layer = 0; layer < features.layers; layer++) {
    for (let i = 0; i < features.squares[layer].length; i++) {
      const square = features.squares[layer][i]
      // Work out the height of the square
      const halfSquareSize = (square.corners.bl.y - square.corners.tl.y)
      if (square.middle.y + halfSquareSize > lowestY) lowestY = square.middle.y
    }
  }
  const bottemYOfPage = 1 / ratio - features.sideBorder
  const subSubSquareSize = features.squareSize / 4
  // How many subSubSquares can we fit in the space between the lowestY and the bottomYOfPage?
  const numberOfSubSubSquares = Math.floor((bottemYOfPage - lowestY) / subSubSquareSize)

  const transformLine = (line, fromSquare, toSquare, scaleFactor) => {
    const translateAndScale = (point) => ({
      x: (point.x - fromSquare.middle.x) / scaleFactor + toSquare.middle.x,
      y: (point.y - fromSquare.middle.y) / scaleFactor + toSquare.middle.y
    })

    return {
      p1: translateAndScale(line.p1),
      p2: translateAndScale(line.p2)
    }
  }

  features.lowestY = lowestY
  features.subSubSquareSize = subSubSquareSize
  features.numberOfSubSubSquaresDown = numberOfSubSubSquares
  features.numberOfSubSubSquaresAcross = Math.floor((1 - (features.sideBorder * 2)) / subSubSquareSize)

  for (let sy = 0; sy < features.numberOfSubSubSquaresDown; sy++) {
    for (let sx = 0; sx < features.numberOfSubSubSquaresAcross; sx++) {
      try {
        const layer = Math.floor(R.prng() * features.layers)
        const squareIndex = Math.floor(R.prng() * features.squares[layer].length)
        const square = features.squares[layer][squareIndex]
        const squareSize = square.corners.tr.x - square.corners.tl.x
        const scaleFactor = squareSize / features.subSubSquareSize

        const middleX = features.sideBorder + (sx * features.subSubSquareSize) + (features.subSubSquareSize / 2)
        const middleY = features.lowestY + (sy * features.subSubSquareSize) + (features.subSubSquareSize / 2)
        const newSquare = { lines: [], middle: { x: middleX, y: middleY } }

        square.lines.forEach(line => {
          newSquare.lines.push(transformLine(line, square, newSquare, scaleFactor))
        })

        features.squares[layer].push(newSquare)
      } catch (er) {}
    }
  }

  // Grab the page height
  let pageHeight = 1 / ratio
  // I want to flip ALL the y values, so I'm going to loop through all the squares
  for (let layer = 0; layer < features.layers; layer++) {
    for (let i = 0; i < features.squares[layer].length; i++) {
      const square = features.squares[layer][i]
      // Loop through all the lines in the square
      for (let j = 0; j < square.lines.length; j++) {
        const line = square.lines[j]
        // Flip the y values
        line.p1.y = pageHeight - line.p1.y
        line.p2.y = pageHeight - line.p2.y
      }
    }
  }

  const strippingRectangles = []
  const numberOfRectangles = Math.floor(R.prng() * 5) + 5 // For 2 to 5 rectangles
  const pageWidth = 1 - (features.sideBorder * 2)
  pageHeight = (1 / ratio) - (features.sideBorder * 2)

  const maxWidth = pageWidth / 2
  const minWidth = pageWidth / 7
  const maxHeight = pageHeight / 2
  const minHeight = pageHeight / 7

  for (let i = 0; i < numberOfRectangles; i++) {
    const width = minWidth + R.prng() * (maxWidth - minWidth)
    const height = minHeight + R.prng() * (maxHeight - minHeight)

    const maxX = pageWidth - width
    const maxY = pageHeight - height
    let x = features.sideBorder + R.prng() * maxX
    let y = features.sideBorder + R.prng() * maxY

    // Move the rectangle to the nearest border 75% of the time
    const distances = [
      { border: 'top', distance: y - features.sideBorder },
      { border: 'bottom', distance: maxY - y },
      { border: 'left', distance: x - features.sideBorder },
      { border: 'right', distance: maxX - x }
    ]

    const nearestBorder = distances.reduce((a, b) => a.distance < b.distance ? a : b).border
    if (R.prng() < 0.75) {
      if (nearestBorder === 'top') y = features.sideBorder
      if (nearestBorder === 'bottom') y = 1 / ratio - features.sideBorder - height
      if (nearestBorder === 'left') x = features.sideBorder
      if (nearestBorder === 'right') x = 1 - features.sideBorder - width
    }

    strippingRectangles.push({ x, y, width, height, nearestBorder })

    // Now we want to remove any lines that have either point inside the rectangle
    // the chance of the line getting removed is based on how far the line is from the edge that's closest to the nearestBorder edge
    // and the opposite edge. So if something is on the right border, lines on the right of the rectange will have a near 100% chance
    // to be removed and the lines on the left will have a near 0% chance to be removed.
    // So first work out which edge is the nearest to the border and which is the opposite edge
    let nearestEdge = null
    let oppositeEdge = null
    let removingDirection = null
    if (nearestBorder === 'top') {
      nearestEdge = y
      oppositeEdge = y + height
      removingDirection = 'down'
    }
    if (nearestBorder === 'bottom') {
      nearestEdge = y + height
      oppositeEdge = y
      removingDirection = 'up'
    }
    if (nearestBorder === 'left') {
      nearestEdge = x
      oppositeEdge = x + width
      removingDirection = 'right'
    }
    if (nearestBorder === 'right') {
      nearestEdge = x + width
      oppositeEdge = x
      removingDirection = 'left'
    }
    // Now loop through all the lines in all the squares finding ones where either point is inside the rectangle, we'll do that first, then
    // work out the chance of it being removed
    for (let layer = 0; layer < features.layers; layer++) {
      for (let i = 0; i < features.squares[layer].length; i++) {
        const square = features.squares[layer][i]
        // Loop through all the lines in the square
        for (let j = 0; j < square.lines.length; j++) {
          const line = square.lines[j]
          // Check if either point is inside the rectangle
          try {
            if ((line.p1.x > x && line.p1.x < x + width && line.p1.y > y && line.p1.y < y + height) || (line.p2.x > x && line.p2.x < x + width && line.p2.y > y && line.p2.y < y + height)) {
              // Grab the midpoint of the line
              const midPoint = {
                x: (line.p1.x + line.p2.x) / 2,
                y: (line.p1.y + line.p2.y) / 2
              }
              // Based on the midPoint work out the chance of it being removed
              let chance = 0
              if (removingDirection === 'down') chance = (midPoint.y - nearestEdge) / (oppositeEdge - nearestEdge)
              if (removingDirection === 'up') chance = (nearestEdge - midPoint.y) / (nearestEdge - oppositeEdge)
              if (removingDirection === 'left') chance = (midPoint.x - nearestEdge) / (oppositeEdge - nearestEdge)
              if (removingDirection === 'right') chance = (nearestEdge - midPoint.x) / (nearestEdge - oppositeEdge)
              // Now remove the line based on the chance
              if (R.prng() > chance) {
                square.lines.splice(j, 1)
                j--
              }
            }
            // If either end of the line is within 0.01 of the edge of the page, we also remove it
            if (line.p1.x < 0.01 || line.p1.x > 0.99 || line.p1.y < 0.01 || line.p1.y > (0.99 / ratio) || line.p2.x < 0.01 || line.p2.x > 0.99 || line.p2.y < 0.01 || line.p2.y > (0.99 / ratio)) {
              square.lines.splice(j, 1)
              j--
            }
          } catch (er) {}
        }
      }
    }
  }

  // Work out the actual Y break point which is the distanceAboveYBreak + the top border
  const yBreakPoint = Math.floor(R.prng() * features.squaresDown) * features.squareSize + features.topBottomBorder
  const xBreakPoint = Math.floor(R.prng() * features.squaresDown) * features.squareSize + features.sideBorder
  features.yBreak = R.prng() < 0.4
  features.xBreak = R.prng() < 0.4
  // const yBreakPoint = 0.5 * (1 / ratio)
  const downShift = (1 / ratio) - yBreakPoint - features.topBottomBorder
  const upShift = yBreakPoint - features.topBottomBorder
  const rightShift = 1 - xBreakPoint - features.sideBorder
  const leftShift = xBreakPoint - features.sideBorder
  for (let layer = 0; layer < features.layers; layer++) {
    for (let i = 0; i < features.squares[layer].length; i++) {
      const square = features.squares[layer][i]
      // Loop through all the lines in the square
      for (let j = 0; j < square.lines.length; j++) {
        const line = square.lines[j]
        // Work out the middle of the line
        line.midPoint = {
          x: (line.p1.x + line.p2.x) / 2,
          y: (line.p1.y + line.p2.y) / 2
        }
        if (features.yBreak) {
          if (line.midPoint.y > yBreakPoint) {
            line.p1.y -= upShift
            line.p2.y -= upShift
          } else {
            line.p1.y += downShift
            line.p2.y += downShift
          }
        }
        if (features.xBreak) {
          if (line.midPoint.x > xBreakPoint) {
            line.p1.x -= leftShift
            line.p2.x -= leftShift
          } else {
            line.p1.x += rightShift
            line.p2.x += rightShift
          }
        }
      }
    }
  }

  // Now delete 50% of the squares
  const squaresToDelete = {}
  let maxSquaresToDelete = 0
  // Loop through the layers and count the squares
  for (let layer = 0; layer < features.layers; layer++) {
    maxSquaresToDelete += features.squares[layer].length
  }
  // Now work out how many to delete
  maxSquaresToDelete = Math.floor(maxSquaresToDelete * 0.3)
  if (features.squaresAcross === 3) maxSquaresToDelete *= 0.333
  if (features.squaresAcross === 5) maxSquaresToDelete *= 0.666

  // Now keep deleting squares until we have deleted enough randomly from random layers
  let squaresToDeleteCount = 0
  while (squaresToDeleteCount < maxSquaresToDelete) {
    // Pick a random layer
    const layer = Math.floor(R.prng() * features.layers)
    // Pick a random square
    const squareIndex = Math.floor(R.prng() * features.squares[layer].length)
    const index = `${layer}__${squareIndex}`
    // Add it to the squaresToDelete array if it's not already there
    if (!squaresToDelete[index]) {
      squaresToDelete[index] = true
      squaresToDeleteCount++
    }
  }
  // Count the number of keys in the squaresToDelete array
  // Empty all the squares in the features object
  for (let layer = 0; layer < features.layers; layer++) {
    const squares = JSON.parse(JSON.stringify(features.squares[layer]))
    features.squares[layer].length = 0
    // Loop through the squares adding them back in, unless they are in the squaresToDelete array
    for (let i = 0; i < squares.length; i++) {
      if (!squaresToDelete[`${layer}__${i}`]) features.squares[layer].push(squares[i])
    }
  }
  // Put the
  // Now loop through the squaresToDelete array and delete the squares
  /*
  for (let i = 0; i < squaresToDelete.length; i++) {
    const squareToDelete = squaresToDelete[i]
    features.squares[squareToDelete.layer].splice(squareToDelete.squareIndex, 1)
  }
  */

  features.strippingRectangles = strippingRectangles

  // Now feed that into the $fx.features object
  const readableFeaturesObj = {}
  readableFeaturesObj['Pallette Name'] = features.palette.name
  readableFeaturesObj.Size = features.squaresAcross
  readableFeaturesObj.layers = features.layers
  readableFeaturesObj.Gridding = 'Normal'
  if (features.cropLines === 0.8) readableFeaturesObj.Gridding = 'Enforced'
  if (features.cropLines === 1) readableFeaturesObj.Gridding = 'Extreme Enforced'
  readableFeaturesObj.Blitting = 'No'
  if (features.yBreak || features.xBreak) readableFeaturesObj.Blitting = 'Yes'
  if (features.yBreak && features.xBreak) readableFeaturesObj.Blitting = 'Extra'

  window.alba.setMetadata(readableFeaturesObj)
  // Drop the features object into the console so we can see it
  console.table(readableFeaturesObj)
}
// Call the setup function straight away, we want this to happen as soon as possible so
// the fxhash system has access to the $fx.features object right away
setup()

/*
 * This is the draw function, it gets called whenever we need to draw the design
 * generall you want to keep all random choices out of here, everything you want
 * to display has already been decided. This function is just drawing what we already
 * know we want to draw.
 *
 * If you need some amount of randomness (for noise, textures, or small details)
 * we'll cover that in a future YouTube video.
 *
 * For the moment we're focusing on a seperation between "data" and "display", here's
 * a few ways of thinking about it...
 * setup() = decisions, drawCanvas() = display
 * setup() = data,      drawCanvas() = markup
 * setup() = backend,   drawCanvas() = frontend
 */
const drawCanvas = async () => {
  // Cancel the next animation frame (we don't really need to do this here, but it's good practice,
  // we don't want to end up having multiple animation frames running at the same time)
  window.cancelAnimationFrame(nextFrame)

  // Grab all the canvas stuff
  const canvas = document.getElementById('target')
  const ctx = canvas.getContext('2d')
  const w = canvas.width
  const h = canvas.height

  /* **************************************************************************
   *
   * This is where your own drawing code goes
   *
   * vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv */

  // Fill in the background in white
  ctx.fillStyle = features.backgroundColour
  ctx.fillRect(0, 0, w, h)

  // Now I want to add texture to the back, first we save the state
  ctx.save()
  // Now we move the origin to the middle of the canvas
  ctx.translate(w / 2, h / 2)
  // Now rotate it 45 degrees
  ctx.rotate(Math.PI / 4)

  // Now we're going to make lots of vertical lines made up of small circles
  ctx.strokeStyle = features.textureColour
  ctx.lineWidth = w / 1000
  //
  let alphaMod = 1.75
  if (features.textureColour === 'white') alphaMod = 9

  // We want to make sure we cover the whole thing so start at -w * 2 to the left and then go to w * 2 to the right
  for (let x = -w * 2; x < w * 2; x += w / 100) {
    // Now we need to move down the page with another loop so we can draw little circles
    // all the way down the page
    for (let y = -h * 2; y < h * 2; y += w / 300) {
      // Now draw the circle
      // the radiusMod can be anything between 400 and 1200
      const radiusMod = 1000 + R.prng() * 2000
      // The alpha can be anything from 0.05 to 0.00
      const alpha = R.prng() * 0.02 * alphaMod
      ctx.globalAlpha = alpha
      ctx.beginPath()
      ctx.arc(x, y, w / radiusMod, 0, Math.PI * 2)
      ctx.stroke()
    }
  }

  // Now are going to draw some concentric circles, made up of lots of little circles, like before
  // Lets have three of them
  for (let i = 0; i < 3; i++) {
    // Pick a random x,y for the middle
    const x = R.prng() * w - (w / 2)
    const y = R.prng() * h - (h / 2)
    // Pick a random radius somewhere between 0.5 and 1.0 of the height
    let radius = h * (0.5 + R.prng() * 0.5) * 0.5
    // Draw a circle filled in with the background colour
    ctx.globalAlpha = 1
    ctx.fillStyle = features.backgroundColour
    ctx.beginPath()
    ctx.arc(x, y, radius, 0, Math.PI * 2)
    ctx.fill()

    // The step is going to be w / 100
    const step = w / 100
    ctx.strokeStyle = features.textureColour
    // Now draw circles until the radius is less than step
    while (radius > step) {
      // Work out the circumference of the circle
      const circumference = Math.PI * radius * 2
      // if a step around the circle is w / 300, how many, rounded, steps are there around the circle?
      const numberOfSteps = Math.round(circumference / (w / 300))
      // Now loop through the steps
      for (let j = 0; j < numberOfSteps; j++) {
        // Work out the angle
        const angle = (Math.PI * 2) / numberOfSteps * j
        // Now work out the x and y for the circle
        const circleX = Math.cos(angle) * radius + x
        const circleY = Math.sin(angle) * radius + y
        // Now draw the circle
        const radiusMod = 1000 + R.prng() * 2000
        // The alpha can be anything from 0.05 to 0.00
        const alpha = R.prng() * 0.02 * alphaMod
        ctx.globalAlpha = alpha
        ctx.beginPath()
        ctx.arc(circleX, circleY, w / radiusMod, 0, Math.PI * 2)
        ctx.stroke()
      }
      // Now reduce the radius by the step amount
      radius -= step
    }
  }
  // Set the alpha back to 1
  ctx.globalAlpha = 1

  // restore the state
  ctx.restore()

  // If the palette is 'Print Me' set the blend mode to multiply
  if (features.palette.name === 'Print Me') ctx.globalCompositeOperation = 'multiply'
  if (features.palette.name === 'Shibuya Reverse') ctx.globalCompositeOperation = 'multiply'
  if (features.palette.name === 'Notaiesque') ctx.globalCompositeOperation = 'multiply'

  if (features.palette.name === 'Shibuya') ctx.globalCompositeOperation = 'screen'
  if (features.palette.name === 'Matrix') ctx.globalCompositeOperation = 'screen'
  if (features.palette.name === 'Oil') ctx.globalCompositeOperation = 'screen'
  if (features.palette.name === 'Ox') ctx.globalCompositeOperation = 'screen'

  // Draw the squares so I can see what they look like
  ctx.lineWidth = w / 10 * (features.squareSize * 0.2)
  // Rounded joins and caps please
  ctx.lineJoin = 'round'
  // ctx.lineCap = 'round'
  for (let layer = 0; layer < features.layers; layer++) {
    for (let i = 0; i < features.squares[layer].length; i++) {
    // Draw the square scaling everything up by the height of the canvas
      const square = features.squares[layer][i]
      // Now draw the lines in the correct colour
      ctx.strokeStyle = features.lineHolder[layer].colour
      ctx.beginPath()
      // Loop through the lines
      for (let j = 0; j < square.lines.length; j++) {
        const line = square.lines[j]
        ctx.moveTo(line.p1.x * w, line.p1.y * w)
        ctx.lineTo(line.p2.x * w, line.p2.y * w)
      }
      ctx.stroke()
    }
  }

  // Set the composite operation back to normal
  ctx.globalCompositeOperation = 'source-over'

  // Write the number of lines in the first layer as text in the middle of the canvas
  /*
  ctx.fillStyle = '#FF0000'
  ctx.font = `${w / 10}px sans-serif`
  ctx.textAlign = 'center'
  ctx.textBaseline = 'middle'
  ctx.fillText(features.lineHolder[0].lines.length, w / 2, h / 2)
  */

  /* ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
   *
   * Above is where your own drawing code goes, everything below here
   * is just handling the animation, and downloading the canvas
   *
   * **************************************************************************/

  // If we haven't taken a thumbnail yet, then take one now
  if (!thumbnailTaken) {
    window.alba.setComplete(true)
    thumbnailTaken = true
  }

  // If we are forcing download, then do that now
  if ('forceDownload' in urlParams && forceDownloaded === false) {
    forceDownloaded = true
    await autoDownloadCanvas()
    // Tell the parent window that we have downloaded, by posting a 'forceDownloaded' message
    // (This is very optional!!)
    window.parent.postMessage('forceDownloaded', '*')
    // Reload the page in 50ms with the forceDownload param still in the url
    // This is so we can download the next one
    setTimeout(() => {
      window.location.href = `${window.location.href}`
    })
  }

  // Draw everything again in the next animation frame, if we are animating
  if (animated) {
    nextFrame = window.requestAnimationFrame(drawCanvas)
  }
}

// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// Below are common functions that don't need to change, they handle starting
// the project, laying out the canvas, and handling downloading snapshots.
//
// You don't need to touch anything below here, see README.md for more info
//
// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

/*
 * This is the init function, it gets called when the DOM is ready. This only ever
 * gets called once to set up the event listeners. Then it kicks things off
 * by calling layoutCanvas(), which in turn calls drawCanvas()
 */
const init = async () => {
  // This is an event listener that gets called when the window is resized. We put it into
  // a timeout so it gets called 100ms _after_ the window has stopped resizing. This is
  // to stop it getting called too often _as_ the window is resized.
  window.addEventListener('resize', async () => {
    //  When we resize we need to layout the canvas again for the new size
    clearTimeout(resizeTmr)
    resizeTmr = setTimeout(async () => {
      await layoutCanvas()
    }, 100)
  })

  // Handle all the keypresses here
  document.addEventListener('keypress', async (e) => {
    e = e || window.event
    // Save the canvas as a PNG
    if (e.key === 's') autoDownloadCanvas()
  })

  //  Now call layout the canvas, which will in turn call drawCanvas()
  await layoutCanvas()
}

/*
 * This function lays out the canvas, and calls drawCanvas() to draw the design
 * This gets called when the window is resized, and when the page first loads.
 *
 * It destroys any existing canvas elements, and creates a new one designed to
 * fit the window size, unless we are forcing the width via the url,
 * in which case it creates a canvas of that width.
 */
const layoutCanvas = async (windowObj = window, urlParamsObj = urlParams) => {
  //  Kill the next animation frame (note, this isn't always used, only if we're animating)
  windowObj.cancelAnimationFrame(nextFrame)

  //  Get the window size, and devicePixelRatio
  const { innerWidth: wWidth, innerHeight: wHeight, devicePixelRatio = 1 } = windowObj
  let dpr = devicePixelRatio
  let cWidth = wWidth
  let cHeight = cWidth / ratio

  // If the height is too big, then we need to adjust the width to fit the height instead
  if (cHeight > wHeight) {
    cHeight = wHeight
    cWidth = wHeight * ratio
  }

  // Grab any canvas elements so we can delete them
  const canvases = document.getElementsByTagName('canvas')
  Array.from(canvases).forEach(canvas => canvas.remove())

  // Now set the target width and height
  let targetHeight = cHeight
  let targetWidth = targetHeight * ratio

  // If we are forcing the width, then use that, and set the dpr to 1
  // (as we want to render at the exact size)
  if ('forceWidth' in urlParams) {
    targetWidth = parseInt(urlParams.forceWidth)
    targetHeight = Math.floor(targetWidth / ratio)
    dpr = 1
  }

  // Update based on the dpr
  targetWidth *= dpr
  targetHeight *= dpr

  // Create a new canvas element, and append it to the body
  // based on all the size stuff we just worked out
  const canvas = document.createElement('canvas')
  canvas.id = 'target'
  canvas.width = Math.floor(targetWidth)
  canvas.height = Math.floor(targetHeight)
  document.body.appendChild(canvas)

  // Now we need to scale the canvas via CSS to make it fit the window
  canvas.style.position = 'absolute'
  canvas.style.width = `${Math.floor(cWidth)}px`
  canvas.style.height = `${Math.floor(cHeight)}px`
  canvas.style.left = `${(wWidth - cWidth) / 2}px`
  canvas.style.top = `${(wHeight - cHeight) / 2}px`

  // Finally we draw the canvas!
  drawCanvas()
}

/*
 * This function converts the canvas to a PNG and downloads it
 * It gets called when the user presses 's', or when we are told to via the URL
 */
const autoDownloadCanvas = async () => {
  const canvas = document.getElementById('target')

  // Create a download link, if we are forcing the id then we add that to the filename
  // i.e. ?forceId=1 will add _0001 to the filename, a url full of "overrides"
  // may look like this:
  // ?forceId=1&forceWidth=2000&forceDownload=true&fxhash=ooABCDEF1234567890
  const element = document.createElement('a')
  const filename = 'forceId' in urlParams
    ? `${prefix}_${urlParams.forceId.toString().padStart(4, '0')}_${window.alba.params.seed}`
    : `${prefix}_${window.alba.params.seed}`
  element.setAttribute('download', filename)

  // Hide the link element
  element.style.display = 'none'
  document.body.appendChild(element)

  // Convert canvas to Blob and set it as the link's href
  const imageBlob = await new Promise(resolve => canvas.toBlob(resolve, 'image/png'))
  element.setAttribute('href', window.URL.createObjectURL(imageBlob))

  // Trigger the download
  element.click()

  // Clean up by removing the link element
  document.body.removeChild(element)
}

/*
 * When everything in the DOM is loaded then we start everything off by calling init()
 *
 * If you have more complicated things going on, like pre-loading images or fonts in
 * some clever way, then you'd do that here, and then call init() when you know
 * everything is ready. For the moment though this one-liner is all we need.
 */
document.addEventListener('DOMContentLoaded', init)
