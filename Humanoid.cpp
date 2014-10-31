#include <Simbody.h>
#include "Humanoid.h"

using namespace std;
using namespace SimTK;

#define ANIMATE

namespace {
const Real D2R = Pi/180; // convert degrees to radians

// Set this <1 for slow motion, >1 for speedup.
const Real RealTimeFactor 
    = 1; // try to run in real time
//  = 0.2; // 5X slower than real time

// Joint stops (coordinate limit forces).
const Real StopStiffness = 1000/D2R; // N-m/deg -> N-m/radian
const Real StopDissipation = 1;

const Real TransitionVelocity = .03; // slide->stick velocity
const Real mu_s = 5;       // Friction coefficients.
const Real mu_d = 5;
const Real mu_v = 0;

// Rubber for feet
//const Real rubber_density = 1100.;  // kg/m^3
//const Real rubber_young   = 0.01e9; // pascals (N/m)
//const Real rubber_poisson = 0.5;    // ratio
//const Real rubber_planestrain = 
//    ContactMaterial::calcPlaneStrainStiffness(rubber_young,rubber_poisson);
//const Real rubber_dissipation = /*0.005*/1;
//
//const ContactMaterial rubber(rubber_planestrain,rubber_dissipation,
//                               mu_s,mu_d,mu_v);
//
//// Concrete for ground
//const Real concrete_density = 2300.;  // kg/m^3
//const Real concrete_young   = 25e9;  // pascals (N/m)
//const Real concrete_poisson = 0.15;    // ratio
//const Real concrete_planestrain = 
//    ContactMaterial::calcPlaneStrainStiffness(concrete_young,concrete_poisson);
//const Real concrete_dissipation = 0.005;
//
//const ContactMaterial concrete(concrete_planestrain,concrete_dissipation,
//                               mu_s,mu_d,mu_v);
//
//const Real ContactSphereRadius = 1.5*.01;
//// Use this clique for contact surfaces on the humanoid that you don't want
//// to have bump into one another. They will still contact Ground.
//const ContactCliqueId ModelClique = ContactSurface::createNewContactClique();

}

Humanoid::Humanoid() 
:   system(), matter(system), forces(system), tracker(system),
    contact(system,tracker), viz(system), userInput(0)
{
//    contact.setTransitionVelocity(TransitionVelocity);
    Force::Gravity(forces,matter,-YAxis,9.80660000);

    // Add the Ground contact geometry. Contact half space has -XAxis normal
    // (right hand wall) so we have to rotate.
//    matter.updGround().updBody()
//        .addContactSurface(Transform(Rotation(-Pi/2,ZAxis),Vec3(0)),
//                    ContactSurface(ContactGeometry::HalfSpace(),concrete));

    constructSystem();
    setUpVisualizer();
}

void Humanoid::findContactStatus(const State& s, 
                                 bool& left, bool& right) const {
    //Real fLeft, fRight;
    //findContactForces(s,fLeft,fRight);
    //left = fLeft > 0; right = fRight > 0;
    left = false;
    right = false;
    for (int i = 1; i < 6; i++)
    {
        left = left || _leftContacts[i]->isEnabled(s);
        right = right || _rightContacts[i]->isEnabled(s);
    }

}

void Humanoid::findContactForces(const State& s, 
                                 Real& fLeft, Real& fRight) const {
    //const int nContacts = contact.getNumContactForces(s);
    //const ContactSnapshot& snapshot = tracker.getActiveContacts(s);
    ////printf("  Humanoid::findContactForces(t=%g) %d contacts active\n", 
    ////    s.getTime(), nContacts);

    //fLeft = fRight = 0;
    //for (int i=0; i < nContacts; ++i) {
    //    const ContactForce& force = contact.getContactForce(s,i);
    //    const ContactId     id    = force.getContactId();
    //    assert(snapshot.hasContact(id));
    //    const Contact&      contact = snapshot.getContactById(id);
    //    const MobilizedBody& b1 = tracker.getMobilizedBody(contact.getSurface1());
    //    const MobilizedBody& b2 = tracker.getMobilizedBody(contact.getSurface2());
    //    const bool left = isLeftFoot(b1) || isLeftFoot(b2);
    //    const bool right = isRightFoot(b1) || isRightFoot(b2);
    //    //printf("    contact %d: bodies %d,%d (left=%d, right=%d)\n", (int)id,
    //    //    (int)b1.getMobilizedBodyIndex(), (int)b2.getMobilizedBodyIndex(),
    //    //    (int)left, (int)right);

    //    if (left)
    //        fLeft += force.getForceOnSurface2()[1].norm();
    //    else if (right)
    //        fRight += force.getForceOnSurface2()[1].norm();
    //}
}


void Humanoid::constructSystem() {
    // Original OpenSim ellipsoid_center.vtp half dimensions.
    const Vec3 ectr(.03, .12, .03); 
    // Original OpenSim sphere.vtp half dimension.
    const Real rad = .5;
    // Original OpenSim block.vtp half dimensions.
    const Vec3 blk(.05,.05,.05);

    // Some convenient rotation matrices for relabeling axes.
    const Rotation Rzmxmy(BodyRotationSequence,  Pi/2, XAxis,  Pi/2, ZAxis);
    const Rotation Rzxy(BodyRotationSequence,   -Pi/2, XAxis, -Pi/2, ZAxis);
    const Rotation Ryzx(BodyRotationSequence,    Pi/2, YAxis,  Pi/2, ZAxis);

    DecorativeSphere orgMarker(0.04*rad);
    orgMarker.setColor(Vec3(.1,.1,.1));


    //--------------------------------------------------------------------------
    //                          Body information 
    //--------------------------------------------------------------------------


    const Real trunkMass = 19.733716299458440;
    const Vec3 trunkCOM(0,0,0);
    Body trunkInfo(MassProperties(trunkMass, trunkCOM,
                                  Inertia(0.355666710204554,
                                          0.224533650416368, 
                                          0.281526683481324, 
                                          0.046737850487895, 
                                          0.000655032101243, 
                                         -0.000572934744554)));

    trunkInfo.addDecoration(Transform(Rotation(.01, ZAxis),
                                      Vec3(-0.025,-0.06,0)),
        DecorativeEllipsoid(ectr.elementwiseMultiply(Vec3(1.6,1.25,1.6)))
            .setColor(White).setOpacity(1));
    trunkInfo.addDecoration(Transform(Rotation(-.3, ZAxis),
                                      Vec3(-0.01,0.16,0)),
        DecorativeEllipsoid(ectr.elementwiseMultiply(Vec3(0.85,0.55,0.85)))
            .setColor(White).setOpacity(1));
    trunkInfo.addDecoration(Transform(Ryzx,Vec3(-0.025,0.09,0)),
        DecorativeEllipsoid(ectr.elementwiseMultiply(Vec3(1,1.43,1)))
            .setColor(White).setOpacity(1));
    trunkInfo.addDecoration(Vec3(0),orgMarker);

    const Real headMass = 4.340524803328146;
    const Vec3 headCOM(0.041488177946752, 0.085892319249215, 0);
    Body headInfo(MassProperties(headMass, headCOM,
                                 Inertia( 0.020576741740382,  
                                          0.016008984554380,  
                                          0.025555859085964,  
                                         -0.004554656543977, 
                                         -0.000103058383929, 
                                          0.000186029116753)
                                 .shiftFromMassCenter(-headCOM,headMass)));
    headInfo.addDecoration(Vec3(0,.11,0),
        DecorativeEllipsoid(Vec3(0.2,0.22,0.2)/2.)
            .setColor(White).setOpacity(1));
    headInfo.addDecoration(Vec3(0),orgMarker);

    const Real pelvisMass = 13.924855817213411;
    const Vec3 pelvisCOM(0.036907589663647, -0.142772863495411, 0);
    Body pelvisInfo(MassProperties(pelvisMass, pelvisCOM,
                                   Inertia( 0.172382614643800,  
                                            0.137961114411544,  
                                            0.128551359933154,  
                                           -0.010239461806632, 
                                            0.001027963710884, 
                                           -0.003286514395970)
                                   .shiftFromMassCenter(-pelvisCOM,pelvisMass)));
    pelvisInfo.addDecoration(Transform(Ryzx,Vec3(0.02,-0.1,0.0)),
        DecorativeEllipsoid(ectr.elementwiseMultiply(Vec3(3.3,1.5,2.4)))
            .setColor(White).setOpacity(1));
    pelvisInfo.addDecoration(Vec3(0),orgMarker);

    // COM z has opposite sign left to right.
    const Real upperarmMass = 2.070989783095760;
    const Vec3 upperarm_rCOM(0.003289136233947,-0.078058926824158,0.065606556342984);
    Body upperarm_rInfo(MassProperties(upperarmMass, upperarm_rCOM,
                      Inertia(0.014082101994221, 
                              0.003604979502976, 
                              0.013769380880709)
                      .shiftFromMassCenter(-upperarm_rCOM, upperarmMass)));
    upperarm_rInfo.addDecoration(Transform(
                    Rotation(BodyRotationSequence,-0.47,XAxis,0.13,ZAxis),
                    Vec3(0.017,-0.12,0.06)),
        DecorativeEllipsoid(ectr.elementwiseMultiply(Vec3(1,1.2,1)))
            .setColor(White).setOpacity(1));
    upperarm_rInfo.addDecoration(Vec3(0),orgMarker);

    const Vec3 upperarm_lCOM(0.003289136233947,-0.078058926824158,-0.065606556342984);
    Body upperarm_lInfo(MassProperties(upperarmMass, upperarm_lCOM,
                      Inertia(0.014082101994221, 
                              0.003604979502976, 
                              0.013769380880709)
                      .shiftFromMassCenter(-upperarm_lCOM, upperarmMass)));
    upperarm_lInfo.addDecoration(Transform(
                    Rotation(BodyRotationSequence,0.47,XAxis,0.13,ZAxis),
                    Vec3(0.017,-0.12,-0.06)),
        DecorativeEllipsoid(ectr.elementwiseMultiply(Vec3(1,1.2,1)))
            .setColor(White).setOpacity(1));
    upperarm_lInfo.addDecoration(Vec3(0),orgMarker);

    // Some signs differ left to right.
    const Real lowerarmMass = 1.106702647320712;
    const Vec3 lowerarm_rCOM(0.031656703591848,-0.089369993258598,0.017231110378866);
    Body lowerarm_rInfo(MassProperties(lowerarmMass, lowerarm_rCOM,
                            Inertia( 0.003846276658463,  
                                     0.001704523106360,  
                                     0.004819186789386,  
                                     0.001553953681336, 
                                    -0.000083971410109, 
                                     0.000083971410109)
                            .shiftFromMassCenter(-lowerarm_rCOM,lowerarmMass)));
    lowerarm_rInfo.addDecoration(Transform(
                    Rotation(BodyRotationSequence,-0.05,XAxis,0.45,ZAxis),
                    Vec3(0.053,-0.111,0.01)),
        DecorativeEllipsoid(ectr.elementwiseMultiply(Vec3(0.95,1.15,0.95)))
            .setColor(White).setOpacity(1));
    lowerarm_rInfo.addDecoration(Vec3(0),orgMarker);

    const Vec3 lowerarm_lCOM(0.031656703591848,-0.089369993258598,-0.017231110378866);
    Body lowerarm_lInfo(MassProperties(lowerarmMass, lowerarm_lCOM,
                            Inertia( 0.003846276658463,  
                                     0.001704523106360,  
                                     0.004819186789386,  
                                     0.001553953681336, 
                                     0.000083971410109, 
                                    -0.000083971410109)
                            .shiftFromMassCenter(-lowerarm_lCOM,lowerarmMass)));
    lowerarm_lInfo.addDecoration(Transform(
                    Rotation(BodyRotationSequence,0.05,XAxis,0.45,ZAxis),
                    Vec3(0.053,-0.111,-0.01)),
        DecorativeEllipsoid(ectr.elementwiseMultiply(Vec3(0.95,1.15,0.95)))
            .setColor(White).setOpacity(1));
    lowerarm_lInfo.addDecoration(Vec3(0),orgMarker);

    // Some signs differ left to right.
    const Real handMass = 0.340742469583528;
    const Vec3 hand_rCOM(0.031681557027587,-0.041582042351409,-0.008872831097566);
    Body hand_rInfo(MassProperties(handMass, hand_rCOM,
                            Inertia( 0.000294382529694,
                                     0.000262531305170, 
                                     0.000326233754218, 
                                     0.000061772071805, 
                                     0.000054050562829, 
                                    -0.000036677167634)
                            .shiftFromMassCenter(-hand_rCOM,handMass)));
    hand_rInfo.addDecoration(Transform(
                    Rotation(0.8,ZAxis),
                    Vec3(0.025,-0.025,0)),
        DecorativeEllipsoid(ectr.elementwiseMultiply(Vec3(0.7,0.3,0.7)))
            .setColor(White).setOpacity(1));
    hand_rInfo.addDecoration(Vec3(0),orgMarker);

    const Vec3 hand_lCOM(0.031681557027587,-0.041582042351409, 0.008872831097566);
    Body hand_lInfo(MassProperties(handMass, hand_lCOM,
                            Inertia( 0.000294382529694,
                                     0.000262531305170, 
                                     0.000326233754218, 
                                     0.000061772071805, 
                                    -0.000054050562829, 
                                     0.000036677167634)
                            .shiftFromMassCenter(-hand_lCOM,handMass)));
    hand_lInfo.addDecoration(Transform(
                    Rotation(0.8,ZAxis),
                    Vec3(0.025,-0.025,0)),
        DecorativeEllipsoid(ectr.elementwiseMultiply(Vec3(0.7,0.3,0.7)))
            .setColor(White).setOpacity(1));
    hand_lInfo.addDecoration(Vec3(0),orgMarker);

    const Real thighMass = 8.082407914884000;
    const Vec3 thigh_rCOM(0,-0.178920728716523,0.001605747837523);
    Body thigh_rInfo(MassProperties(thighMass, thigh_rCOM,
                            Inertia( 0.116351777130644,
                                     0.030499980412887, 
                                     0.122695077900275 )
                            .shiftFromMassCenter(-thigh_rCOM,thighMass)));
    thigh_rInfo.addDecoration(Transform(
                    Rotation(BodyRotationSequence, -.01,XAxis, -.012,ZAxis),
                    Vec3(-0.002,-0.205,0.003)),
        DecorativeEllipsoid(ectr.elementwiseMultiply(Vec3(1.3,1.9,1.3)))
            .setColor(White).setOpacity(1));
    thigh_rInfo.addDecoration(Vec3(0),orgMarker);

    const Vec3 thigh_lCOM(0,-0.178920728716523,-0.001605747837523);
    Body thigh_lInfo(MassProperties(thighMass, thigh_lCOM,
                            Inertia( 0.116351777130644,
                                     0.030499980412887, 
                                     0.122695077900275 )
                            .shiftFromMassCenter(-thigh_lCOM,thighMass)));
    thigh_lInfo.addDecoration(Transform(
                    Rotation(BodyRotationSequence, .01,XAxis, -.012,ZAxis),
                    Vec3(-0.002,-0.205,-0.003)),
        DecorativeEllipsoid(ectr.elementwiseMultiply(Vec3(1.3,1.9,1.3)))
            .setColor(White).setOpacity(1));
    thigh_lInfo.addDecoration(Vec3(0),orgMarker);


    const Real shankMass = 3.222323418816000;
    const Vec3 shank_rCOM(0,-0.182765070363067,0.005552190835500);
    Body shank_rInfo(MassProperties(shankMass, shank_rCOM,
                            Inertia( 0.043804477493817,
                                     0.004432595936874, 
                                     0.044412873014564 )
                            .shiftFromMassCenter(-shank_rCOM,shankMass)));
    shank_rInfo.addDecoration(Transform(
                    Rotation(0.030,XAxis),
                    Vec3(0.0,-0.21,-0.005)),
        DecorativeEllipsoid(ectr.elementwiseMultiply(Vec3(1.27,1.8,1.27)))
            .setColor(White).setOpacity(1));
    shank_rInfo.addDecoration(Vec3(0),orgMarker);

    const Vec3 shank_lCOM(0,-0.182765070363067,-0.005552190835500);
    Body shank_lInfo(MassProperties(shankMass, shank_lCOM,
                            Inertia( 0.043804477493817,
                                     0.004432595936874, 
                                     0.044412873014564 )
                            .shiftFromMassCenter(-shank_lCOM,shankMass)));
    shank_lInfo.addDecoration(Transform(
                    Rotation(-0.030,XAxis),
                    Vec3(0.0,-0.21,0.005)),
        DecorativeEllipsoid(ectr.elementwiseMultiply(Vec3(1.27,1.8,1.27)))
            .setColor(White).setOpacity(1));
    shank_lInfo.addDecoration(Vec3(0),orgMarker);

    const Real footMass = 1.172905458264000;
    const Vec3 foot_rCOM(0.035606945567853,-0.051617802456029,-0.000574057583573);
    const Inertia foot_Ic(0.001313654113256, // central inertia
                          0.003659465029784, 
                          0.003847129903106);
    const Real Ifac = 1; // for playing with foot inertia; should be 1
    Body foot_rInfo(MassProperties(footMass, foot_rCOM, 
                      (Ifac*foot_Ic).shiftFromMassCenter(-foot_rCOM,footMass)));
    foot_rInfo.addDecoration(Vec3(0.052,-0.043,0.0),
        DecorativeBrick(blk.elementwiseMultiply(Vec3(1.65,0.6,0.8)))
            .setColor(White).setOpacity(1));
    foot_rInfo.addDecoration(Vec3(0),orgMarker);

    const Vec3 foot_lCOM(0.035606945567853,-0.051617802456029,0.000574057583573);
    Body foot_lInfo(MassProperties(footMass, foot_lCOM,
                      (Ifac*foot_Ic).shiftFromMassCenter(-foot_lCOM,footMass)));
    foot_lInfo.addDecoration(Vec3(0.052,-0.043,0.0),
        DecorativeBrick(blk.elementwiseMultiply(Vec3(1.65,0.6,0.8)))
            .setColor(White).setOpacity(1));
    foot_lInfo.addDecoration(Vec3(0),orgMarker);


    const Real toesMass = 0.20;
    // This inertia looks artificially large.
    const Inertia toes_Ic(0.087132432150,0.0174264864299,0.0871324321496);
    const Vec3 toes_rCOM(0.023716003435794,-0.001184334594970,-0.002484544347746);
    Body toes_rInfo(MassProperties(toesMass, toes_rCOM,
                            toes_Ic.shiftFromMassCenter(-toes_rCOM,toesMass)));
    toes_rInfo.addDecoration(Vec3(0.03,0.014,0.0),
        DecorativeBrick(blk.elementwiseMultiply(Vec3(0.6,0.3,0.8)))
            .setColor(White).setOpacity(1));
    toes_rInfo.addDecoration(Vec3(0),orgMarker);

    const Vec3 toes_lCOM(0.023716003435794,-0.001184334594970,0.002484544347746);
    Body toes_lInfo(MassProperties(toesMass, toes_lCOM,
                            toes_Ic.shiftFromMassCenter(-toes_lCOM,toesMass)));
    toes_lInfo.addDecoration(Vec3(0.03,0.014,0.0),
        DecorativeBrick(blk.elementwiseMultiply(Vec3(0.6,0.3,0.8)))
            .setColor(White).setOpacity(1));
    toes_lInfo.addDecoration(Vec3(0),orgMarker);

    /*
    // Add contact spheres to feet and toes.
    ContactSurface contactBall(ContactGeometry::Sphere(ContactSphereRadius), 
                               rubber);
    contactBall.joinClique(ModelClique);
    DecorativeSphere contactArt(ContactSphereRadius);
    contactArt.setColor(Magenta);

    const Real footsphereht = -.07;
    for (int x=0; x <= 1; ++x)
        for (int z=0; z <= 1; ++z) {
            const Vec3 rctr(heelball[x], footsphereht,  rlatmed[z]);
            const Vec3 lctr(heelball[x], footsphereht, -rlatmed[z]);
            foot_rInfo.addContactSurface(rctr, contactBall);
            foot_rInfo.addDecoration(rctr, contactArt);
            foot_lInfo.addContactSurface(lctr, contactBall);
            foot_lInfo.addDecoration(lctr, contactArt);
        }

    // Balls just at toe tips.
    const Real toetip = .06;
    for (int z=0; z <= 1; ++z) {
        const Vec3 rctr(toetip, 0,  rlatmed[z]);
        const Vec3 lctr(toetip, 0, -rlatmed[z]);
        toes_rInfo.addContactSurface(rctr, contactBall);
        toes_rInfo.addDecoration(rctr, contactArt);
        toes_lInfo.addContactSurface(lctr, contactBall);
        toes_lInfo.addDecoration(lctr, contactArt);
    }
    */


    //--------------------------------------------------------------------------
    //                         Mobilized Bodies
    //--------------------------------------------------------------------------
    trunk = MobilizedBody::Free(matter.updGround(), Vec3(0),
                                trunkInfo,          Vec3(0));

    // Neck angles are: q0=extension (about z), q1=bending (x), q2=rotation (y).
    head = MobilizedBody::ORIENT(
        trunk,    Transform(Rzxy, Vec3(0.010143822053248,0.222711680750785,0)),
        headInfo, Rzxy);
    Force::MobilityLinearStop(forces, head, MobilizerQIndex(0), //extension
        StopStiffness,StopDissipation, -80*D2R, 50*D2R);
    Force::MobilityLinearStop(forces, head, MobilizerQIndex(1), //bending
        StopStiffness,StopDissipation, -60*D2R, 50*D2R);
    Force::MobilityLinearStop(forces, head, MobilizerQIndex(2), //rotation
        StopStiffness,StopDissipation, -80*D2R, 80*D2R);

    // Back angles are: q0=tilt (about z), q1=list (x), q2=rotation (y).
    pelvis = MobilizedBody::ORIENT(
        trunk, Transform(Rzxy, Vec3(-0.019360589663647,-0.220484136504589,0)),
        pelvisInfo, Rzxy);
    Force::MobilityLinearStop(forces, pelvis, MobilizerQIndex(0), //tilt
        StopStiffness,StopDissipation, -5*D2R, 10*D2R);
    Force::MobilityLinearStop(forces, pelvis, MobilizerQIndex(1), //list
        StopStiffness,StopDissipation, -5*D2R, 5*D2R);
    Force::MobilityLinearStop(forces, pelvis, MobilizerQIndex(2), //rotation
        StopStiffness,StopDissipation, -15*D2R, 15*D2R);

    // Right arm.
    //-----------
    // Shoulder angles are q0=flexion, q1=adduction, q2=rotation
    upperarm_r = MobilizedBody::ORIENT(
        trunk, Transform(Rzxy, 
                Vec3(-0.023921136233947,0.079313926824158,0.164710443657016)),
        upperarm_rInfo, Rzxy);
    Force::MobilityLinearStop(forces, upperarm_r, MobilizerQIndex(0), //flexion
        StopStiffness,StopDissipation, -80*D2R, 160*D2R);
    Force::MobilityLinearStop(forces, upperarm_r, MobilizerQIndex(1), //adduction
        StopStiffness,StopDissipation, -45*D2R, 45*D2R); // TODO: was -90:90
    Force::MobilityLinearStop(forces, upperarm_r, MobilizerQIndex(2), //rotation
        StopStiffness,StopDissipation, -20*D2R, 20*D2R);

    // Elbow angles are q0=flexion, q1=rotation
    lowerarm_r = MobilizedBody::Universal(
        upperarm_r, Transform(Rotation(-Pi/2,YAxis),
                Vec3(0.033488432642100,-0.239093933565560,0.117718445964118)),
        lowerarm_rInfo, Rotation(-Pi/2,YAxis));
    Force::MobilityLinearStop(forces, lowerarm_r, MobilizerQIndex(0), //flexion
        StopStiffness,StopDissipation, 0*D2R, 120*D2R);
    Force::MobilityLinearStop(forces, lowerarm_r, MobilizerQIndex(1), //rotation
        StopStiffness,StopDissipation, -90*D2R, 40*D2R);

    hand_r = MobilizedBody::Weld(
        lowerarm_r, Vec3(0.110610146564261,-0.232157950907188,0.014613941476432),
        hand_rInfo, Vec3(0));

    // Left arm.
    //----------
    upperarm_l = MobilizedBody::ORIENT(
        trunk, Transform(Rzmxmy, 
                Vec3(-0.023921136233947,0.079313926824158,-0.164710443657016)),
        upperarm_lInfo, Rzmxmy);
    Force::MobilityLinearStop(forces, upperarm_l, MobilizerQIndex(0), //flexion
        StopStiffness,StopDissipation, -80*D2R, 160*D2R);
    Force::MobilityLinearStop(forces, upperarm_l, MobilizerQIndex(1), //adduction
        StopStiffness,StopDissipation, -45*D2R, 45*D2R); // TODO: was -90:90
    Force::MobilityLinearStop(forces, upperarm_l, MobilizerQIndex(2), //rotation
        StopStiffness,StopDissipation, -20*D2R, 20*D2R);

    lowerarm_l = MobilizedBody::Universal(
        upperarm_l, Transform(
                Rotation(BodyRotationSequence, Pi/2,YAxis, Pi,ZAxis),
                Vec3(0.033488432642100,-0.239093933565560,-0.117718445964118)),
        lowerarm_lInfo, Rotation(BodyRotationSequence, Pi/2,YAxis, Pi,ZAxis));
    Force::MobilityLinearStop(forces, lowerarm_l, MobilizerQIndex(0), //flexion
        StopStiffness,StopDissipation, 0*D2R, 120*D2R);
    Force::MobilityLinearStop(forces, lowerarm_l, MobilizerQIndex(1), //rotation
        StopStiffness,StopDissipation, -90*D2R, 40*D2R);

    hand_l = MobilizedBody::Weld(
        lowerarm_l, Vec3(0.110610146564261,-0.232157950907188,-0.014613941476432),
        hand_lInfo, Vec3(0));


    // Right leg.
    //-----------
    // Hip angles are q0=flexion, q1=adduction, q2=rotation.
    thigh_r = MobilizedBody::ORIENT(
        pelvis, Transform(Rzxy, 
                Vec3(0.029343644095793,-0.180413783750097,0.117592252162477)),
        thigh_rInfo, Rzxy);
    Force::MobilityLinearStop(forces, thigh_r, MobilizerQIndex(0), //flexion
        StopStiffness,StopDissipation, -60*D2R, 165*D2R);
    Force::MobilityLinearStop(forces, thigh_r, MobilizerQIndex(1), //adduction
        StopStiffness,StopDissipation, -20*D2R, 20*D2R);
    Force::MobilityLinearStop(forces, thigh_r, MobilizerQIndex(2), //rotation
        StopStiffness,StopDissipation, -120*D2R, 20*D2R);

    // Knee angle is q0=extension
    shank_r = MobilizedBody::Pin(
        thigh_r, Vec3(-0.005,-0.416780050422019,0.004172557002023),
        shank_rInfo, Vec3(0));
    Force::MobilityLinearStop(forces, shank_r, MobilizerQIndex(0), //extension
        StopStiffness,StopDissipation, -165*D2R, 0*D2R);

    // Ankle angles are q0=dorsiflexion, q1=inversion.
    foot_r = MobilizedBody::Universal(
        shank_r, Transform(
            Rotation(BodyRotationSequence,Pi/2,XAxis,Pi,YAxis,Pi/2,ZAxis),
            Vec3(0,-0.420937226867266,-0.011971751580927)),
        foot_rInfo, 
            Rotation(BodyRotationSequence,Pi/2,XAxis,Pi,YAxis,Pi/2,ZAxis));
    Force::MobilityLinearStop(forces, foot_r, MobilizerQIndex(0), //dorsiflexion
        StopStiffness,StopDissipation, -50*D2R, 30*D2R);
    Force::MobilityLinearStop(forces, foot_r, MobilizerQIndex(1), //inversion
        StopStiffness,StopDissipation, -2*D2R, 35*D2R);

    // Toe angle is q0=dorsiflexion
    toes_r = MobilizedBody::Pin(
        foot_r, Vec3(0.134331942132059,-0.071956467861059,-0.000000513235827),
        toes_rInfo, Vec3(0));
    Force::MobilityLinearStop(forces, toes_r, MobilizerQIndex(0), //dorsiflexion
        StopStiffness,StopDissipation, 0*D2R, 30*D2R);

    // Left leg.
    //----------
    thigh_l = MobilizedBody::ORIENT(
        pelvis, Transform(Rzmxmy, 
                Vec3(0.029343644095793,-0.180413783750097,-0.117592252162477)),
        thigh_lInfo, Rzmxmy);
    Force::MobilityLinearStop(forces, thigh_l, MobilizerQIndex(0), //flexion
        StopStiffness,StopDissipation, -60*D2R, 165*D2R);
    Force::MobilityLinearStop(forces, thigh_l, MobilizerQIndex(1), //adduction
        StopStiffness,StopDissipation, -20*D2R, 20*D2R);
    Force::MobilityLinearStop(forces, thigh_l, MobilizerQIndex(2), //rotation
        StopStiffness,StopDissipation, -120*D2R, 20*D2R);

    shank_l = MobilizedBody::Pin(
        thigh_l, Vec3(-0.005,-0.416780050422019,-0.004172557002023),
        shank_lInfo, Vec3(0));
    Force::MobilityLinearStop(forces, shank_l, MobilizerQIndex(0), //extension
        StopStiffness,StopDissipation, -165*D2R, 0*D2R);

    foot_l = MobilizedBody::Universal(
        shank_l, Transform(
            Rotation(BodyRotationSequence,-Pi/2,XAxis,Pi,YAxis,-Pi/2,ZAxis),
            Vec3(0,-0.420937226867266,0.011971751580927)),
        foot_lInfo, 
            Rotation(BodyRotationSequence,-Pi/2,XAxis,Pi,YAxis,-Pi/2,ZAxis));
    Force::MobilityLinearStop(forces, foot_l, MobilizerQIndex(0), //dorsiflexion
        StopStiffness,StopDissipation, -50*D2R, 30*D2R);
    Force::MobilityLinearStop(forces, foot_l, MobilizerQIndex(1), //inversion
        StopStiffness,StopDissipation, -2*D2R, 35*D2R);

    toes_l = MobilizedBody::Pin (
        foot_l, Vec3(0.134331942132059,-0.071956467861059,0.000000513235827),
        toes_lInfo, Vec3(0));
    Force::MobilityLinearStop(forces, toes_l, MobilizerQIndex(0), //dorsiflexion
        StopStiffness,StopDissipation, 0*D2R, 30*D2R);

    //--------------------------------------------------------------------------
    //                         Rigid contact
    //--------------------------------------------------------------------------
//    matter.updGround().updBody().addContactSurface(Transform(NegX2Y, Vec3(0)),
//            ContactSurface(ContactGeometry::HalfSpace(), ContactMaterial(ContactMaterial::calcPlaneStrainStiffness(2.5e9, 0.4), 0, 0, 0, 0));
    const Real CoefRest = 0;
    const Real mu_s = 0.8;
    const Real mu_d = 0.5;
    const Real mu_v = 0.0;
    const Vec2 rlatmed(.04,-.04), heelball(-.02,.125);
    const Real footsphereht = -.07;
    _leftContacts.resize(6);
    _rightContacts.resize(6);
    int i = 0;
    for (int x=0; x<=1; ++x) {
        for (int z=0; z<=1; ++z) {
            const Vec3 rctr(heelball[x], footsphereht,  rlatmed[z]);
            const Vec3 lctr(heelball[x], footsphereht, -rlatmed[z]);
            PointPlaneContact * ppcr = new PointPlaneContact(
                    matter.updGround(), YAxis, 0.,
                    foot_r, rctr, CoefRest, mu_s, mu_d, mu_v);
            PointPlaneContact * ppcl = new PointPlaneContact(
                    matter.updGround(), YAxis, 0.,
                    foot_l, lctr, CoefRest, mu_s, mu_d, mu_v);
            matter.adoptUnilateralContact(ppcr);
            matter.adoptUnilateralContact(ppcl);

            _rightContacts[i] = ppcr;
            _leftContacts[i] = ppcl;

            i += 1;
        }
    }
    i = 0;

    const Real toetip = .06;
    for (int z=0; z <= 1; ++z) {
        const Vec3 rctr(toetip, 0,  rlatmed[z]);
        const Vec3 lctr(toetip, 0, -rlatmed[z]);
        PointPlaneContact * ppcr = new PointPlaneContact(
                matter.updGround(), YAxis, 0.,
                toes_r, rctr, CoefRest, mu_s, mu_d, mu_v);
        PointPlaneContact * ppcl = new PointPlaneContact(
                matter.updGround(), YAxis, 0.,
                toes_l, lctr, CoefRest, mu_s, mu_d, mu_v);
        matter.adoptUnilateralContact(ppcr);
        matter.adoptUnilateralContact(ppcl);

        _rightContacts[i] = ppcr;
        _leftContacts[i] = ppcl;

        i += 1;
    }



}

void Humanoid::setUpVisualizer() {
    // Set up visualization and ask for a frame every 1/30 second.
    viz.setShowSimTime(true); viz.setShowFrameRate(true);
    userInput = new Visualizer::InputSilo();
    viz.addInputListener(userInput);   
    #ifdef ANIMATE
    system.addEventReporter(new Visualizer::Reporter(viz, RealTimeFactor/30));
    #endif
    DecorativeText help("Any input to start; ESC to quit");
    help.setIsScreenText(true);
    viz.addDecoration(MobilizedBodyIndex(0),Vec3(0),help);
    matter.setShowDefaultGeometry(false);
}

void Humanoid::fillInActuatorMap(const State& s) {
    act2coord.resize(NumActuators); 

    act2coord[neck_extension]= getQUIx(s,head, 0); 
    act2coord[neck_bending]  = getQUIx(s,head, 1); 
    act2coord[neck_rotation] = getQUIx(s,head, 2); 

    act2coord[back_tilt]     = getQUIx(s,pelvis, 0); 
    act2coord[back_list]     = getQUIx(s,pelvis, 1); 
    act2coord[back_rotation] = getQUIx(s,pelvis, 2); 

    act2coord[shoulder_r_flexion]   = getQUIx(s,upperarm_r, 0); 
    act2coord[shoulder_r_adduction] = getQUIx(s,upperarm_r, 1); 
    act2coord[shoulder_r_rotation]  = getQUIx(s,upperarm_r, 2);

    act2coord[elbow_r_flexion]  = getQUIx(s,lowerarm_r, 0); 
    act2coord[elbow_r_rotation] = getQUIx(s,lowerarm_r, 1); 

    act2coord[shoulder_l_flexion]   = getQUIx(s,upperarm_l, 0); 
    act2coord[shoulder_l_adduction] = getQUIx(s,upperarm_l, 1); 
    act2coord[shoulder_l_rotation]  = getQUIx(s,upperarm_l, 2);

    act2coord[elbow_l_flexion]  = getQUIx(s,lowerarm_l, 0); 
    act2coord[elbow_l_rotation] = getQUIx(s,lowerarm_l, 1); 

    act2coord[hip_r_flexion]    = getQUIx(s,thigh_r, 0); 
    act2coord[hip_r_adduction]  = getQUIx(s,thigh_r, 1); 
    act2coord[hip_r_rotation]   = getQUIx(s,thigh_r, 2); 

    act2coord[knee_r_extension] = getQUIx(s,shank_r, 0); 

    act2coord[ankle_r_dorsiflexion] = getQUIx(s,foot_r, 0); 
    act2coord[ankle_r_inversion]    = getQUIx(s,foot_r, 1); 

    act2coord[mtp_r_dorsiflexion] = getQUIx(s,toes_r, 0); 

    act2coord[hip_l_flexion]    = getQUIx(s,thigh_l, 0); 
    act2coord[hip_l_adduction]  = getQUIx(s,thigh_l, 1); 
    act2coord[hip_l_rotation]   = getQUIx(s,thigh_l, 2); 

    act2coord[knee_l_extension] = getQUIx(s,shank_l, 0); 

    act2coord[ankle_l_dorsiflexion] = getQUIx(s,foot_l, 0); 
    act2coord[ankle_l_inversion]    = getQUIx(s,foot_l, 1); 

    act2coord[mtp_l_dorsiflexion] = getQUIx(s,toes_l, 0); 
}
