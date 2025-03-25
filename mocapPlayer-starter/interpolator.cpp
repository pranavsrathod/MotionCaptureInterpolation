#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>
#include "motion.h"
#include "interpolator.h"
#include "transform.h"
#include "types.h"
#include <fstream>
#include "performanceCounter.h"
#include <iostream>

Interpolator::Interpolator()
{
  //Set default interpolation type
  m_InterpolationType = LINEAR;

  //set default angle representation to use for interpolation
  m_AngleRepresentation = EULER;
}

Interpolator::~Interpolator()
{
}

//Create interpolated motion
void Interpolator::Interpolate(Motion * pInputMotion, Motion ** pOutputMotion, int N) 
{
  //Allocate new motion
  *pOutputMotion = new Motion(pInputMotion->GetNumFrames(), pInputMotion->GetSkeleton()); 

  //Perform the interpolation
  if ((m_InterpolationType == LINEAR) && (m_AngleRepresentation == EULER))
    LinearInterpolationEuler(pInputMotion, *pOutputMotion, N);
  else if ((m_InterpolationType == LINEAR) && (m_AngleRepresentation == QUATERNION))
    LinearInterpolationQuaternion(pInputMotion, *pOutputMotion, N);
  else if ((m_InterpolationType == BEZIER) && (m_AngleRepresentation == EULER))
    BezierInterpolationEuler(pInputMotion, *pOutputMotion, N);
  else if ((m_InterpolationType == BEZIER) && (m_AngleRepresentation == QUATERNION))
    BezierInterpolationQuaternion(pInputMotion, *pOutputMotion, N);
  else
  {
    printf("Error: unknown interpolation / angle representation type.\n");
    exit(1);
  }
}

void Interpolator::LinearInterpolationEuler(Motion * pInputMotion, Motion * pOutputMotion, int N)
{
  int inputLength = pInputMotion->GetNumFrames(); // frames are indexed 0, ..., inputLength-1

  int startKeyframe = 0;
  count_performance.StartCounter();
  while (startKeyframe + N + 1 < inputLength)
  {
    int endKeyframe = startKeyframe + N + 1;

    Posture * startPosture = pInputMotion->GetPosture(startKeyframe);
    Posture * endPosture = pInputMotion->GetPosture(endKeyframe);

    // copy start and end keyframe
    pOutputMotion->SetPosture(startKeyframe, *startPosture);
    pOutputMotion->SetPosture(endKeyframe, *endPosture);

    // Extract lfemur X-axis rotation (first value in bone_rotation[lfemur])
    // double startAngle = startPosture->bone_rotation[2].p[0];  // lfemur X
    // double endAngle = endPosture->bone_rotation[2].p[0];      // lfemur X

    // interpolate in between
    for(int frame=1; frame<=N; frame++)
    {
      Posture interpolatedPosture;
      double t = 1.0 * frame / (N+1);

      // interpolate root position
      interpolatedPosture.root_pos = startPosture->root_pos * (1-t) + endPosture->root_pos * t;

      // interpolate bone rotations
      for (int bone = 0; bone < MAX_BONES_IN_ASF_FILE; bone++)
        interpolatedPosture.bone_rotation[bone] = startPosture->bone_rotation[bone] * (1-t) + endPosture->bone_rotation[bone] * t;

      pOutputMotion->SetPosture(startKeyframe + frame, interpolatedPosture);
    }

    startKeyframe = endKeyframe;
  }
  count_performance.StopCounter();
  printf("Time taken to perform Linear Euler Interpolation: %f\n", count_performance.GetElapsedTime());
  for(int frame=startKeyframe+1; frame<inputLength; frame++)
    pOutputMotion->SetPosture(frame, *(pInputMotion->GetPosture(frame)));
  
}

void Interpolator::Rotation2Euler(double R[9], double angles[3])
{
  double cy = sqrt(R[0]*R[0] + R[3]*R[3]);

  if (cy > 16*DBL_EPSILON) 
  {
    angles[0] = atan2(R[7], R[8]);
    angles[1] = atan2(-R[6], cy);
    angles[2] = atan2(R[3], R[0]);
  } 
  else 
  {
    angles[0] = atan2(-R[5], R[4]);
    angles[1] = atan2(-R[6], cy);
    angles[2] = 0;
  }

  for(int i=0; i<3; i++)
    angles[i] *= 180 / M_PI;
}

void Interpolator::Euler2Rotation(double angles[3], double R[9])
{
  double rX[4][4], rY[4][4], rZ[4][4], multZY[4][4], xmultZY[4][4];
  rotationX(rX, angles[0]); rotationY(rY, angles[1]); rotationZ(rZ, angles[2]);
  matrix_mult(rZ, rY, multZY); matrix_mult(multZY, rX, xmultZY);
  for(int i=0; i<3; i++){
    for(int j=0; j<3; j++){
      R[j + 3 * i] = xmultZY[i][j];
    }
  }

}


void Interpolator::BezierInterpolationEuler(Motion * pInputMotion, Motion * pOutputMotion, int N)
{
    int inputLength = pInputMotion->GetNumFrames();
    // pOutputMotion->SetNumFrames(inputLength);  // Ensure correct frame storage

    int startKeyframe = 0;
    count_performance.StartCounter();
    while (startKeyframe + N + 1 < inputLength)
    {
        int endKeyframe = startKeyframe + N + 1;

        // Get keyframe postures
        Posture* startPosture = pInputMotion->GetPosture(startKeyframe);
        Posture* endPosture = pInputMotion->GetPosture(endKeyframe);


        // Compute previous and next keyframes for control point calculation
        Posture* prevPosture = (startKeyframe >= N + 1) ? pInputMotion->GetPosture(startKeyframe - (N + 1)) : startPosture;
        Posture* nextPosture = (endKeyframe + N + 1 < inputLength) ? pInputMotion->GetPosture(endKeyframe + (N + 1)) : endPosture;

        // Store keyframes
        pOutputMotion->SetPosture(startKeyframe, *startPosture);
        pOutputMotion->SetPosture(endKeyframe, *endPosture);
        // printf("Keyframe stored: %d and %d\n", startKeyframe, endKeyframe);

        // === Bezier Interpolation Loop ===
        for (int frame = 1; frame <= N; frame++)
        {
            double t = (double)frame / (N + 1);
            Posture interpolatedPosture;  
            // interpolatedPosture.root_pos = startPosture->root_pos * (1 - t) + endPosture->root_pos * t;
            vector p1 = startPosture->root_pos + (endPosture->root_pos - prevPosture->root_pos) * (1.0 / 3.0);
            vector p2 = endPosture->root_pos - (nextPosture->root_pos - startPosture->root_pos) * (1.0 / 3.0);
            interpolatedPosture.root_pos = DeCasteljauEuler(t, startPosture->root_pos, p1, p2, endPosture->root_pos);
            for (int bone = 0; bone < MAX_BONES_IN_ASF_FILE; bone++)
            {
                // Extract Euler angles
                vector anglesStart = startPosture->bone_rotation[bone];
                vector anglesEnd = endPosture->bone_rotation[bone];
                vector anglesPrev = prevPosture->bone_rotation[bone];
                vector anglesNext = nextPosture->bone_rotation[bone];

                // Compute Bezier control points
                vector p1 = anglesStart + (anglesEnd - anglesPrev) * (1.0 / 3.0);
                vector p2 = anglesEnd - (anglesNext - anglesStart) * (1.0 / 3.0);

                // Apply Bezier interpolation using De Casteljau’s algorithm
                vector resultAngles = DeCasteljauEuler(t, anglesStart, p1, p2, anglesEnd);

                // Store the interpolated angles
                interpolatedPosture.bone_rotation[bone] = resultAngles;
            }

            // Store interpolated posture
            pOutputMotion->SetPosture(startKeyframe + frame, interpolatedPosture);
            // printf("Interpolated frame stored: %d\n", startKeyframe + frame);
        }

        startKeyframe = endKeyframe;
    }

    count_performance.StopCounter();
    printf("Time taken to perform Bezier Euler Interpolation: %f\n", count_performance.GetElapsedTime());

    // Copy remaining frames
    for (int frame = startKeyframe + 1; frame < inputLength; frame++)
    {
        pOutputMotion->SetPosture(frame, *(pInputMotion->GetPosture(frame)));
        // printf("Remaining frame stored: %d\n", frame);
    }

    // printf("Total frames in output motion: %d\n", pOutputMotion->GetNumFrames());
}

void Interpolator::LinearInterpolationQuaternion(Motion * pInputMotion, Motion * pOutputMotion, int N)
{
  int inputLength = pInputMotion->GetNumFrames();
  int startKeyframe = 0;
  
  count_performance.StartCounter();
  while (startKeyframe + N + 1 < inputLength)
  {
    int endKeyframe = startKeyframe + N + 1;

    Posture * startPosture = pInputMotion->GetPosture(startKeyframe);
    Posture * endPosture = pInputMotion->GetPosture(endKeyframe);

    // Copy start and end keyframe
    pOutputMotion->SetPosture(startKeyframe, *startPosture);
    pOutputMotion->SetPosture(endKeyframe, *endPosture);

    // Convert root rotation Euler angles to quaternions
    Quaternion<double> qRootStart, qRootEnd, qRootInterpolated;
    Euler2Quaternion(startPosture->bone_rotation[0].p, qRootStart);
    Euler2Quaternion(endPosture->bone_rotation[0].p, qRootEnd);

    // Interpolate in between
    for(int frame = 1; frame <= N; frame++)
    {
      Posture interpolatedPosture;
      double t = 1.0 * frame / (N+1);

      // Interpolate root position (Linear)
      interpolatedPosture.root_pos = startPosture->root_pos * (1 - t) + endPosture->root_pos * t;
      qRootInterpolated = Slerp(t, qRootStart, qRootEnd);
      Quaternion2Euler(qRootInterpolated, interpolatedPosture.bone_rotation[0].p);

      // Interpolate bone rotations using SLERP
      for (int bone = 0; bone < MAX_BONES_IN_ASF_FILE; bone++)
      {
        Quaternion<double> qStart, qEnd, qInterpolated;

        // Ensure `bone_rotation[bone].p` is correctly assigned
        // if (!startPosture->bone_rotation[bone].p || !endPosture->bone_rotation[bone].p)
        // {
        //     printf("ERROR: Bone rotation data is NULL for bone %d\n", bone);
        //     continue;
        // }

        // Convert Euler angles to Quaternions
        Euler2Quaternion(startPosture->bone_rotation[bone].p, qStart);
        Euler2Quaternion(endPosture->bone_rotation[bone].p, qEnd);

        // Debugging print
        // printf("Bone %d: SLERP between QStart(%f, %f, %f, %f) and QEnd(%f, %f, %f, %f)\n", 
        //        bone, qStart.Gets(), qStart.Getx(), qStart.Gety(), qStart.Getz(), 
        //        qEnd.Gets(), qEnd.Getx(), qEnd.Gety(), qEnd.Getz());

        // Perform SLERP
        qInterpolated = Slerp(t, qStart, qEnd);

        // Convert back to Euler angles
        Quaternion2Euler(qInterpolated, interpolatedPosture.bone_rotation[bone].p);
      }

      // Store interpolated posture
      pOutputMotion->SetPosture(startKeyframe + frame, interpolatedPosture);
    }

    startKeyframe = endKeyframe;
  }
  count_performance.StopCounter();
  printf("Time taken to perform SLERP Quaternion Interpolation: %f\n", count_performance.GetElapsedTime());

  for(int frame = startKeyframe + 1; frame < inputLength; frame++)
    pOutputMotion->SetPosture(frame, *(pInputMotion->GetPosture(frame)));
}

void Interpolator::BezierInterpolationQuaternion(Motion * pInputMotion, Motion * pOutputMotion, int N)
{
  int inputLength = pInputMotion->GetNumFrames();
  int startKeyframe = 0;

  while (startKeyframe + N + 1 < inputLength)
  {
      int endKeyframe = startKeyframe + N + 1;
      int prevKeyframe = startKeyframe - (N + 1);
      int nextKeyframe = endKeyframe + (N + 1);

      Posture * startPosture = pInputMotion->GetPosture(startKeyframe);
      Posture * endPosture = pInputMotion->GetPosture(endKeyframe);
      Posture * prevPosture = (prevKeyframe >= 0) ? pInputMotion->GetPosture(prevKeyframe) : startPosture;
      Posture * nextPosture = (nextKeyframe < inputLength) ? pInputMotion->GetPosture(nextKeyframe) : endPosture;

      pOutputMotion->SetPosture(startKeyframe, *startPosture);
      pOutputMotion->SetPosture(endKeyframe, *endPosture);

      count_performance.StartCounter();
      for (int frame = 1; frame <= N; frame++)
      {
          double t = (double)frame / (N + 1);
          Posture interpolatedPosture;
          //interpolatedPosture.root_pos = startPosture->root_pos * (1 - t) + endPosture->root_pos * t;

          vector p1 = startPosture->root_pos + (endPosture->root_pos - prevPosture->root_pos) * (1.0 / 3.0);
          vector p2 = endPosture->root_pos - (nextPosture->root_pos - startPosture->root_pos) * (1.0 / 3.0);
          interpolatedPosture.root_pos = DeCasteljauEuler(t, startPosture->root_pos, p1, p2, endPosture->root_pos);

          Quaternion<double> qStart, qEnd, qPrev, qNext, a_n, b_n, qInterpolated;

          // Convert root rotation Euler angles to quaternions
          Euler2Quaternion(prevPosture->bone_rotation[0].p, qPrev);
          Euler2Quaternion(startPosture->bone_rotation[0].p, qStart);
          Euler2Quaternion(endPosture->bone_rotation[0].p, qEnd);
          Euler2Quaternion(nextPosture->bone_rotation[0].p, qNext);

          // Compute Bezier control points in quaternion space
          Quaternion<double> temp1 = Double(qPrev, qStart);
          Quaternion<double> a_n_bar = Slerp(0.5, temp1, qEnd);
          a_n = Slerp(1.0 / 3.0, qStart, a_n_bar);
          b_n = Slerp(1.0 / 3.0, qEnd, a_n_bar);

          // Perform Bezier SLERP interpolation
          qInterpolated = DeCasteljauQuaternion(t, qStart, a_n, b_n, qEnd);

          // Convert interpolated quaternion back to Euler angles
          Quaternion2Euler(qInterpolated, interpolatedPosture.bone_rotation[0].p);
          
          for (int bone = 0; bone < MAX_BONES_IN_ASF_FILE; bone++)
          {
              Quaternion<double> qPrev, qStart, qEnd, qNext, a_n, b_n, a_n_bar;

              Euler2Quaternion(prevPosture->bone_rotation[bone].p, qPrev);
              Euler2Quaternion(startPosture->bone_rotation[bone].p, qStart);
              Euler2Quaternion(endPosture->bone_rotation[bone].p, qEnd);
              Euler2Quaternion(nextPosture->bone_rotation[bone].p, qNext);

              if (startKeyframe == 0) {
                  // Compute a_1 using formula from Slide 33
                  a_n_bar = Slerp(2.0, qNext, qStart);
                  a_n = Slerp(1.0 / 3.0, qStart, a_n_bar);

                  // Compute b_n normally
                  b_n = Slerp(1.0 / 3.0, qEnd, a_n_bar);
              } 
              else if (endKeyframe == inputLength - 1) {
                  // Compute b_N using formula from Slide 33
                  a_n_bar = Slerp(2.0, qPrev, qEnd);
                  b_n = Slerp(1.0 / 3.0, qEnd, a_n_bar);

                  // Compute a_n normally
                  a_n = Slerp(1.0 / 3.0, qStart, a_n_bar);
              } 
              else {
                  // Compute Bezier control points normally
                  Quaternion<double> temp = Double(qPrev, qStart);
                  a_n_bar = Slerp(0.5, temp, qEnd);
                  a_n = Slerp(1.0 / 3.0, qStart, a_n_bar);
                  b_n = Slerp(1.0 / 3.0, qEnd, a_n_bar);
              }

              Quaternion<double> qInterpolated = DeCasteljauQuaternion(t, qStart, a_n, b_n, qEnd);
              Quaternion2Euler(qInterpolated, interpolatedPosture.bone_rotation[bone].p);
          }

          pOutputMotion->SetPosture(startKeyframe + frame, interpolatedPosture);
      }

      startKeyframe = endKeyframe;
  }
  count_performance.StopCounter();
  printf("Time taken to perform Bezier SLERP Quaternion Interpolation: %f\n", count_performance.GetElapsedTime());

  for (int frame = startKeyframe + 1; frame < inputLength; frame++)
      pOutputMotion->SetPosture(frame, *(pInputMotion->GetPosture(frame)));
}

void Interpolator::Euler2Quaternion(double angles[3], Quaternion<double> & q) 
{
  // students should implement this
  // Convert Euler angles from degrees to radians
    // double radAngles[3] = {
    //     angles[0] * M_PI / 180.0,
    //     angles[1] * M_PI / 180.0,
    //     angles[2] * M_PI / 180.0
    // };

    // Convert Euler angles to a rotation matrix
    // double R[9];  // 3x3 rotation matrix
    // Euler2Rotation(angles, R);  // 🔹 Ensure this function exists

    // // Convert rotation matrix to quaternion
    // q = Quaternion<double>::Matrix2Quaternion(R);
    // q.Normalize();

    // printf("Euler2Quaternion Input: (%f, %f, %f)\n", angles[0], angles[1], angles[2]);

    double R[9];  // 3x3 rotation matrix
    Euler2Rotation(angles, R);  // Ensure this function exists

    q = Quaternion<double>::Matrix2Quaternion(R);
    q.Normalize();

    // printf("Euler2Quaternion Output: (%f, %f, %f, %f)\n", q.Gets(), q.Getx(), q.Gety(), q.Getz());
}

void Interpolator::Quaternion2Euler(Quaternion<double> & q, double angles[3]) 
{
  // students should implement this
  // Convert quaternion to rotation matrix
    double R[9];
    q.Quaternion2Matrix(R);

    // Convert rotation matrix to Euler angles (XYZ order)
    Rotation2Euler(R, angles);  // Ensure this function exists in interpolator.cpp

    // // Convert back to degrees
    // angles[0] *= 180.0 / M_PI;
    // angles[1] *= 180.0 / M_PI;
    // angles[2] *= 180.0 / M_PI;
}

Quaternion<double> Interpolator::Slerp(double t, Quaternion<double> & qStart, Quaternion<double> & qEnd_)
{
    Quaternion<double> result;
    double angle, costheta, sintheta;

    // Compute the dot product (cosTheta)
    costheta = qStart.Gets() * qEnd_.Gets() + 
               qStart.Getx() * qEnd_.Getx() + 
               qStart.Gety() * qEnd_.Gety() + 
               qStart.Getz() * qEnd_.Getz();

    // if (costheta > 0.9995)
    // {
    //     Quaternion<double> result1 = (qStart * (1 - t) + qEnd_ * t);
    //     result1.Normalize();
    //     return result1;
    // }
    //Ensure shortest path in quaternion space
    if (costheta < 0)
    {
        qEnd_ = (-1.0) * qEnd_;  // Flip quaternion
        costheta = -costheta;
    }

    // Compute the angle
    costheta = std::max(-1.0, std::min(1.0, costheta));  // Ensure valid range
    angle = acos(costheta);

    // If the angle is 0, return the start quaternion (no rotation needed)
    if (angle == 0.0)
        return qStart;

    // Compute SLERP using spherical interpolation formula
    sintheta = sin(angle);
    result = (sin((1 - t) * angle) / sintheta) * qStart + 
             (sin(t * angle) / sintheta) * qEnd_;

    return result;
}

Quaternion<double> Interpolator::Double(Quaternion<double> p, Quaternion<double> q)
{
  // students should implement this
  // Quaternion<double> result;
  // return result;
  // return Slerp(2.0, p, q); well This is also valid.
  
  // In Double(p, q) = 2(p*q)q – p, p*q is a dot product. - Tips in aassignment description

  Quaternion<double> result;
  double costheta;
  costheta = p.Gets() * q.Gets() + p.Getx() * q.Getx() + p.Gety() * q.Gety() + p.Getz() * q.Getz();
  result = ((2 * costheta) * q) - p;
  result.Normalize();
  return result;
}

vector Interpolator::DeCasteljauEuler(double t, vector p0, vector p1, vector p2, vector p3)
{
  // students should implement this
  // vector result;
  // return result;
  vector p01 = p0 * (1 - t) + p1 * t;
  vector p12 = p1 * (1 - t) + p2 * t;
  vector p23 = p2 * (1 - t) + p3 * t;

  vector p012 = p01 * (1 - t) + p12 * t;
  vector p123 = p12 * (1 - t) + p23 * t;

  return p012 * (1 - t) + p123 * t;
}

Quaternion<double> Interpolator::DeCasteljauQuaternion(double t, Quaternion<double> p0, Quaternion<double> p1, Quaternion<double> p2, Quaternion<double> p3)
{
  // students should implement this
  Quaternion<double> Q0 = Slerp(t, p0, p1), Q1 = Slerp(t, p1, p2), Q2 = Slerp(t, p2, p3);
  Quaternion<double> R0 = Slerp(t, Q0, Q1), R1 = Slerp(t, Q1, Q2);
  
  Quaternion<double> result = Slerp(t, R0, R1);
  return result;
}

