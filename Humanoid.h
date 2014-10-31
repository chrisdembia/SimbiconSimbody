#ifndef SIMBICON_SIMBODY_HUMANOID_H_
#define SIMBICON_SIMBODY_HUMANOID_H_

#include <Simbody.h>
#include <vector>

// This model should use Gimbal joints because it depends on qdot=u. But you
// can test what happens with Ball joints but picking them here.
//#define USE_BALLS
#ifdef USE_BALLS
    #define ORIENT Ball
#else
    #define ORIENT Gimbal
#endif

//==============================================================================
// Build the SIMBICON 3D Humanoid model
//==============================================================================

// These are the actuators in the order they were defined in the original
// model. These are indices into the controls Vector; the controls must then
// be mapped into the generalized forces Vector which is not ordered the same.
enum Actuator {
    neck_extension = 0,
    neck_bending = 1,
    neck_rotation = 2,
    back_tilt = 3,
    back_list = 4,
    back_rotation = 5,
    shoulder_r_flexion = 6,
    shoulder_r_adduction = 7,
    shoulder_r_rotation = 8,
    elbow_r_flexion = 9,
    elbow_r_rotation = 10,
    shoulder_l_flexion = 11,
    shoulder_l_adduction = 12,
    shoulder_l_rotation = 13,
    elbow_l_flexion = 14,
    elbow_l_rotation = 15,
    hip_r_adduction = 16,
    hip_r_flexion = 17,
    hip_r_rotation = 18,
    knee_r_extension = 19,
    ankle_r_inversion = 20,
    ankle_r_dorsiflexion = 21,
    mtp_r_dorsiflexion = 22,
    hip_l_adduction = 23,
    hip_l_flexion = 24,
    hip_l_rotation = 25,
    knee_l_extension = 26,
    ankle_l_inversion = 27,
    ankle_l_dorsiflexion = 28,
    mtp_l_dorsiflexion = 29,

    NumActuators = mtp_l_dorsiflexion + 1
};
class Humanoid {
public:
    Humanoid();

    // You must call this before using the model.
    void fillInActuatorMap(const SimTK::State& s);

    SimTK::QIndex getQIndex(Actuator act) const
    {   assert(!act2coord.empty()); return act2coord[act].first; }
    SimTK::UIndex getUIndex(Actuator act) const
    {   assert(!act2coord.empty()); return act2coord[act].second; }

    bool isRightFoot(const SimTK::MobilizedBody& mobod) const 
    {   return     mobod.isSameMobilizedBody(foot_r) 
                || mobod.isSameMobilizedBody(toes_r); }
    bool isLeftFoot(const SimTK::MobilizedBody& mobod) const 
    {   return     mobod.isSameMobilizedBody(foot_l) 
                || mobod.isSameMobilizedBody(toes_l); }

    // Return the number of active (force producing) contacts.
    int getNumContacts(const SimTK::State& s) const
    {   return contact.getNumContactForces(s); }

    // Set left,right to true if the corresponding foot is in contact.
    void findContactStatus(const SimTK::State& s, bool& left, bool& right) const;

    // Set fLeft, fRight to total force magnitude due to contact on that foot, 
    // zero if no contact.
    void findContactForces(const SimTK::State& s, 
                           SimTK::Real& fLeft, SimTK::Real& fRight) const;

    SimTK::MultibodySystem              system;
    SimTK::SimbodyMatterSubsystem       matter;
    SimTK::GeneralForceSubsystem        forces;
    SimTK::ContactTrackerSubsystem      tracker;
    SimTK::CompliantContactSubsystem    contact;
    SimTK::Visualizer                   viz;
    SimTK::Visualizer::InputSilo*       userInput; // just a ref; not owned here

    // Mobilized bodies.
    SimTK::MobilizedBody::Free          trunk;
    SimTK::MobilizedBody::ORIENT        head;
    SimTK::MobilizedBody::ORIENT        pelvis;
    SimTK::MobilizedBody::ORIENT        upperarm_r;
    SimTK::MobilizedBody::Universal     lowerarm_r;
    SimTK::MobilizedBody::Weld          hand_r;
    SimTK::MobilizedBody::ORIENT        upperarm_l;
    SimTK::MobilizedBody::Universal     lowerarm_l;
    SimTK::MobilizedBody::Weld          hand_l;
    SimTK::MobilizedBody::ORIENT        thigh_r;
    SimTK::MobilizedBody::Pin           shank_r;
    SimTK::MobilizedBody::Universal     foot_r;
    SimTK::MobilizedBody::Pin           toes_r;
    SimTK::MobilizedBody::ORIENT        thigh_l;
    SimTK::MobilizedBody::Pin           shank_l;
    SimTK::MobilizedBody::Universal     foot_l;
    SimTK::MobilizedBody::Pin           toes_l;
private:
    void constructSystem();
    void setUpVisualizer();

    static std::pair<SimTK::QIndex,SimTK::UIndex> getQUIx
       (const SimTK::State& s, const SimTK::MobilizedBody& mobod,int which)
    {   SimTK::QIndex q0=mobod.getFirstQIndex(s); 
        SimTK::UIndex u0=mobod.getFirstUIndex(s); 
        return std::make_pair(SimTK::QIndex(q0+which), SimTK::UIndex(u0+which)); 
    }

    // Use this to find indices to coordinate q and speed u associated with
    // each actuator.
    SimTK::Array_<std::pair<SimTK::QIndex,SimTK::UIndex>, int> act2coord;
    std::vector<SimTK::PointPlaneContact*> _leftContacts;
    std::vector<SimTK::PointPlaneContact*> _rightContacts;
};



#endif // SIMBICON_SIMBODY_HUMANOID_H_
