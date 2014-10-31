#include "RigidContact.h"

//#undef NDEBUG

using namespace SimTK; 
//==============================================================================
//                        IMPACT HANDLING (CONTACT ON)
//==============================================================================

//------------------------------ HANDLE EVENT ----------------------------------
// There are three different witness functions that cause this handler to get
// invoked. The position- and velocity-level ones require an impulse. The
// acceleration-level one just requires recalculating the active set, in the
// same manner as liftoff or friction transition events.

void ContactOn::
handleEvent(State& s, Real accuracy, bool& shouldTerminate) const 
{
    SimTK_DEBUG3("\nhandle %d impact@%.17g (%s)\n", m_which, s.getTime(),
         m_stage.getName().c_str());

    if (m_stage == Stage::Acceleration) {
        SimTK_DEBUG("\n---------------CONTACT ON (ACCEL)--------------\n");
        SimTK_DEBUG2("CONTACT triggered by element %d @t=%.15g\n", 
            m_which, s.getTime());
        m_mbs.realize(s, Stage::Acceleration);

        #ifndef NDEBUG
        std::cout << " triggers=" << s.getEventTriggers() << "\n";
        #endif

        m_unis.selectActiveConstraints(s, accuracy);
        SimTK_DEBUG("---------------CONTACT ON (ACCEL) done.--------------\n");
        return;
    }

    MyElementSubset proximal;
    m_unis.findProximalElements(s, accuracy, proximal);

    // Zero out accumulated impulses and perform any other necessary 
    // initialization of contact and friction elements.
    for (unsigned i=0; i < proximal.m_contact.size(); ++i) {
        const int which = proximal.m_contact[i];
        m_unis.updContactElement(which)
            .initializeForImpact(s, m_unis.getCaptureVelocity());
    }
    for (unsigned i=0; i < proximal.m_friction.size(); ++i) {
        const int which = proximal.m_friction[i];
        m_unis.updFrictionElement(which).initializeForStiction(s);
    }

    SimTK_DEBUG("\n---------------------CONTACT ON---------------------\n");
    SimTK_DEBUG3("\nIMPACT (%s) for conta2t %d at t=%.16g\n", 
        m_stage.getName().c_str(), m_which, s.getTime());
    SimTK_DEBUG2("  %d/%d proximal contact/friction elements\n", 
        proximal.m_contact.size(), proximal.m_friction.size());

    m_unis.showConstraintStatus(s, "ENTER IMPACT (CONTACT ON)");

    bool needMoreCompression = true;
    while (needMoreCompression) {
        processCompressionPhase(proximal, s);
        needMoreCompression = false;
        if (processExpansionPhase(proximal, s)) {
            for (unsigned i=0; i < proximal.m_contact.size(); ++i) {
                const int which = proximal.m_contact[i];
                const MyContactElement& contact = 
                    m_unis.getContactElement(which);
                if (contact.getVerr(s) < 0) {
                    needMoreCompression = true;
                    break;
                }
            }
        }
    }

    // Record new previous slip velocities for all the sliding friction
    // since velocities have changed. First loop collects the velocities.
    m_mbs.realize(s, Stage::Velocity);
    for (unsigned i=0; i < proximal.m_friction.size(); ++i) {
        const int which = proximal.m_friction[i];
        MyFrictionElement& fric = m_unis.updFrictionElement(which);
        if (!fric.isMasterActive(s) || fric.isSticking(s)) continue;
        fric.recordSlipDir(s);
    }

    // Now update all the previous slip direction state variables from the
    // recorded values.
    for (unsigned i=0; i < proximal.m_friction.size(); ++i) {
        const int which = proximal.m_friction[i];
        const MyFrictionElement& fric = m_unis.getFrictionElement(which);
        if (!fric.isMasterActive(s) || fric.isSticking(s)) continue;
        fric.updatePreviousSlipDirFromRecorded(s);
    }

    m_unis.selectActiveConstraints(s, accuracy);

    SimTK_DEBUG("\n-----------------END CONTACT ON---------------------\n");
}




//------------------------ PROCESS COMPRESSION PHASE ---------------------------
// Given a list of proximal unilateral constraints (contact and stiction),
// determine which ones are active in the least squares solution for the
// constraint impulses. Candidates are those constraints that meet the 
// kinematic proximity condition -- for contacts, position less than
// tolerance; for stiction, master contact is proximal. Also, any
// constraint that is currently active is a candidate, regardless of its
// kinematics.
//
// TODO: NOT ENABLING STICTION AT ALL; assuming stiction on initially leads
// to a lot of wasted time and solution difficulty. Must start with sliding.
// TODO: sliding friction impulses
//
// Algorithm
// ---------
// We're given a set of candidates including contacts and stiction. If any are
// inactive, activate them.
// -- at this point all verr==0, some impulses f might be < 0
//
// loop
// - Calculate impulses with the current active set
//     (iterate sliding impulses until f=mu_d*N to tol, where N is normal 
//      impulse)
// - Calculate f for active constraints, verr for inactive
// - If all f>=0, verr>=0 -> break loop
// - Check for verr < 0 [tol?]. Shouldn't happen but if it does must turn on the
//     associated constraint for the worst violation, then -> continue loop
// - Find worst (most negative) offender:
//    contact offense  = fc < 0 ? fc : 0
//    stiction offense = mu_d*max(0, fc) - |fs|
// - Choose constraint to deactivate:
//     worst is a stiction constraint: choose it
//     worst is a contact constraint: if it has stiction, chose that
//                                    otherwise choose the contact constraint
// - Inactivate chosen constraint
//     (if stiction, record impending slip direction & N for stick->slide)
// end loop 
//
void ContactOn::
processCompressionPhase(MyElementSubset&    proximal,
                        State&              s) const
{
    const int ReviseNormalNIters = 6;

    SimTK_DEBUG2("Entering processCompressionPhase(): "
        "%d/%d impact/stick candidates ...\n", proximal.m_contact.size(),
        proximal.m_friction.size());

    if (proximal.m_contact.empty()) {
        // Can't be any friction either, if there are no contacts.
        SimTK_DEBUG("EXIT processCompressionPhase: no proximal candidates.\n");
        return;
    }

    Vector lambda;
    const Vector u0 = s.getU(); // save presenting velocity

    // Assume at first that all proximal contacts will participate. This is 
    // necessary to ensure that we get a least squares solution for the impulse 
    // involving as many constraints as possible sharing the impulse. 
    // TODO: note stiction is unconditionally disabled
    m_unis.enableConstraintSubset(proximal, false, s); 

    int pass=0, nContactsDisabled=0, nStictionDisabled=0, nContactsReenabled=0;
    while (true) {
        ++pass; 
        SimTK_DEBUG1("processCompressionPhase(): pass %d\n", pass);

        // Given an active set, evaluate impulse multipliers & forces, and
        // evaluate resulting constraint velocity errors.
        calcStoppingImpulse(proximal, s, lambda);
        // TODO: ignoring sliding impacts; should converge here.
        updateVelocities(u0, lambda, s);

        // Scan all proximal contacts to find the active one that has the
        // most negative normal force, and the inactive one that has the 
        // most negative velocity error (hopefully none will).

        int worstActiveContact=-1; Real mostNegativeContactImpulse=0;
        int worstInactiveContact=-1; Real mostNegativeVerr=0;
        
        SimTK_DEBUG("analyzing impact constraints ...\n");
        for (unsigned i=0; i < proximal.m_contact.size(); ++i) {
            const int which = proximal.m_contact[i];
            SimTK_DEBUG1("  %d: ", which);
            const MyContactElement& cont = m_unis.getContactElement(which);
            if (cont.isDisabled(s)) {
                const Real verr = cont.getVerr(s);
                SimTK_DEBUG1("off verr=%g\n", verr);
                if (verr < mostNegativeVerr) {
                    worstInactiveContact = which;
                    mostNegativeVerr = verr;
                }
            } else {
                const Real f = 
                    cont.getMyValueFromConstraintSpaceVector(s, lambda);
                SimTK_DEBUG1("on impulse=%g\n", f);
                if (f < mostNegativeContactImpulse) {
                    worstActiveContact = which;
                    mostNegativeContactImpulse = f;
                }
            }
        }

        // This is bad and might cause cycling.
        if (mostNegativeVerr < 0) {
            SimTK_DEBUG2("  !!! Inactive contact %d violated, verr=%g\n", 
                worstInactiveContact, mostNegativeVerr);
            const MyContactElement& cont = 
                m_unis.getContactElement(worstInactiveContact);
            //TODO -- must use a tolerance or prevent looping
            //++nContactsReenabled;
            //cont.enable(s);
            //continue;
        }

        SimTK_DEBUG("  All inactive contacts are satisfied.\n");

        #ifndef NDEBUG
        if (mostNegativeContactImpulse == 0)
            printf("  All active contacts are satisfied.\n");
        else 
            printf("  Active contact %d was worst violator with f=%g\n",
                worstActiveContact, mostNegativeContactImpulse);
        #endif

        int worstActiveStiction=-1; Real mostNegativeStictionCapacity=0;     
        SimTK_DEBUG("analyzing stiction constraints ...\n");
        for (unsigned i=0; i < proximal.m_friction.size(); ++i) {
            const int which = proximal.m_friction[i];
            SimTK_DEBUG1("  %d: ", which);
            const MyFrictionElement& fric = m_unis.getFrictionElement(which);
            if (!fric.isSticking(s)) {
                SimTK_DEBUG("off\n");
                continue;
            }
            // TODO: Kludge -- must be point contact.
            const MyPointContactFriction& ptfric = 
                dynamic_cast<const MyPointContactFriction&>(fric);
            const MyPointContact& cont = ptfric.getMyPointContact();
            const Real N = cont.getMyValueFromConstraintSpaceVector(s, lambda);
            const Real fsmag = 
                ptfric.getMyImpulseMagnitudeFromConstraintSpaceVector(s, lambda);
            const Real mu_d = fric.getDynamicFrictionCoef();
            const Real capacity = mu_d*std::max(N,Real(0)) - fsmag;
            SimTK_DEBUG2("on capacity=%g (N=%g)\n", capacity, N);

            if (capacity < mostNegativeStictionCapacity) {
                worstActiveStiction = which;
                mostNegativeStictionCapacity = capacity;
            }
        }

        #ifndef NDEBUG
        if (mostNegativeStictionCapacity == 0)
            printf("  All active stiction constraints are satisfied.\n");
        else 
            printf("  Active stiction %d was worst violator with capacity=%g\n",
                worstActiveStiction, mostNegativeStictionCapacity);
        #endif

        if (mostNegativeContactImpulse==0 && mostNegativeStictionCapacity==0) {
            SimTK_DEBUG("DONE. Current active set is a winner.\n");
            break;
        }

        // Restore original velocity.
        s.updU() = u0;

        if (mostNegativeStictionCapacity <= mostNegativeContactImpulse) {
            SimTK_DEBUG1("  Disable stiction %d\n", worstActiveStiction);
            MyFrictionElement& fric = 
                m_unis.updFrictionElement(worstActiveStiction);
            // TODO: need the impulse version of this
            //fric.recordImpendingSlipInfo(s);
            ++nStictionDisabled;
            fric.disableStiction(s);
            continue;
        }

        // A contact constraint was the worst violator. If that contact
        // constraint has an active stiction constraint, we have to disable
        // the stiction constraint first.
        SimTK_DEBUG1("  Contact %d was the worst violator.\n", 
            worstActiveContact);
        const MyContactElement& cont = 
            m_unis.getContactElement(worstActiveContact);
        assert(!cont.isDisabled(s));

        if (cont.hasFrictionElement()) {
            MyFrictionElement& fric = cont.updFrictionElement();
            if (fric.isSticking(s)) {
                SimTK_DEBUG1("  ... but must disable stiction %d first.\n",
                    fric.getFrictionIndex());
                // TODO: need the impulse version of this
                //fric.recordImpendingSlipInfo(s);
                ++nStictionDisabled;
                fric.disableStiction(s);
                continue;
            }
        }

        SimTK_DEBUG1("  Disable contact %d\n", worstActiveContact); 
        ++nContactsDisabled;
        cont.disable(s);
        // Go back for another pass.
    }


    // Now update the entries for each proximal constraint to reflect the
    // compression impulse and post-compression velocity.
    SimTK_DEBUG("  Compression results:\n");
    for (unsigned i=0; i < proximal.m_contact.size(); ++i) {
        const int which = proximal.m_contact[i];
        MyContactElement& cont = m_unis.updContactElement(which);
        if (!cont.isDisabled(s))
            cont.recordImpulse(MyContactElement::Compression, s, lambda);
        SimTK_DEBUG4("  %d %3s: Ic=%g, V=%g\n",
            which, cont.isDisabled(s) ? "off" : "ON", 
            cont.getCompressionImpulse(), cont.getVerr(s));
    }

    SimTK_DEBUG("... compression phase done.\n");
}

//------------------------- PROCESS EXPANSION PHASE ----------------------------
bool ContactOn::
processExpansionPhase(MyElementSubset&  proximal,
                      State&            s) const
{
    SimTK_DEBUG("Entering processExpansionPhase() ...\n");

    // Generate an expansion impulse if there were any active contacts that
    // still have some restitution remaining.
    Vector expansionImpulse;

    bool anyChange = false;
    for (unsigned i=0; i<proximal.m_contact.size(); ++i) {
        const int which = proximal.m_contact[i];
        MyContactElement& uni = m_unis.updContactElement(which);
        if (uni.isDisabled(s)||uni.isRestitutionDone()
            ||uni.getEffectiveCoefRest()==0
            ||uni.getCompressionImpulse()<=0)
            continue;
        uni.setMyExpansionImpulse(s, uni.getEffectiveCoefRest(), 
                                  expansionImpulse);
        uni.recordImpulse(MyContactElement::Expansion,s,expansionImpulse);
        uni.setRestitutionDone(true);
        anyChange = true;
    }

    if (!anyChange) {
        SimTK_DEBUG("... no expansion impulse -- done.\n");
        return false;
    }

    // We generated an expansion impulse. Apply it and update velocities.
    updateVelocities(Vector(), expansionImpulse, s);

    // Release any constraint that now has a positive velocity.
    Array_<int> toDisable;
    for (unsigned i=0; i < proximal.m_contact.size(); ++i) {
        const int which = proximal.m_contact[i];
        const MyContactElement& uni = m_unis.getContactElement(which);
        if (!uni.isDisabled(s) && uni.getVerr(s) > 0)
            toDisable.push_back(which);
    }

    // Now do the actual disabling (can't mix this with checking velocities)
    // because disabling invalidates Instance stage.
    for (unsigned i=0; i < toDisable.size(); ++i) {
        const int which = toDisable[i];
        const MyContactElement& uni = m_unis.getContactElement(which);
        uni.disable(s);
    }

    SimTK_DEBUG("  Expansion results:\n");
    m_mbs.realize(s, Stage::Velocity);
    for (unsigned i=0; i < proximal.m_contact.size(); ++i) {
        const int which = proximal.m_contact[i];
        const MyContactElement& uni = m_unis.getContactElement(which);
        SimTK_DEBUG4("  %d %3s: Ie=%g, V=%g\n",
            which, uni.isDisabled(s) ? "off" : "ON", 
            uni.getExpansionImpulse(), uni.getVerr(s));
    }

    SimTK_DEBUG("... expansion phase done.\n");

    return true;
}




//-------------------------- CALC STOPPING IMPULSE -----------------------------
// Calculate the impulse that eliminates all residual velocity for the
// current set of enabled constraints.
// Note: if you have applied impulses also (like sliding friction), 
// convert to generalized impulse f, then to desired delta V in constraint
// space like this: deltaV = G*M\f; add that to the verrs to get the total
// velocity change that must be produced by the impulse.
void ContactOn::
calcStoppingImpulse(const MyElementSubset&  proximal,
                    const State&            s,
                    Vector&                 lambda0) const
{
    const SimbodyMatterSubsystem& matter = m_mbs.getMatterSubsystem();
    m_mbs.realize(s, Stage::Dynamics); // TODO: should only need Position
    Vector desiredDeltaV;  // in constraint space
    SimTK_DEBUG("  Entering calcStoppingImpulse() ...\n");
    bool gotOne = false;
    for (unsigned i=0; i < proximal.m_contact.size(); ++i) {
        const int which = proximal.m_contact[i];
        const MyContactElement& uni = m_unis.getContactElement(which);
        if (uni.isDisabled(s))
            continue;
        SimTK_DEBUG2("    uni constraint %d enabled, v=%g\n",
            which, uni.getVerr(s));
        uni.setMyDesiredDeltaV(s, desiredDeltaV);
        gotOne = true;
    }
    for (unsigned i=0; i < proximal.m_friction.size(); ++i) {
        const int which = proximal.m_friction[i];
        const MyFrictionElement& fric = m_unis.getFrictionElement(which);
        if (!fric.isSticking(s))
            continue;
        SimTK_DEBUG2("    friction constraint %d enabled, |v|=%g\n",
            which, fric.getActualSlipSpeed(s));
        fric.setMyDesiredDeltaV(s, desiredDeltaV);
        gotOne = true;
    }

    if (gotOne) matter.solveForConstraintImpulses(s, desiredDeltaV, lambda0);
    else lambda0.clear();
#ifndef NDEBUG
    std::cout << "  ... done. Stopping impulse=" << lambda0 << std::endl;
#endif
}



//---------------------------- UPDATE VELOCITIES -------------------------------
void ContactOn::
updateVelocities(const Vector& u0, const Vector& lambda, State& state) const {
    const SimbodyMatterSubsystem& matter = m_mbs.getMatterSubsystem();
    Vector f, deltaU;
    assert(u0.size()==0 || u0.size() == state.getNU());

    m_mbs.realize(state, Stage::Dynamics); // TODO: Position
    matter.multiplyByGTranspose(state,lambda,f);
    matter.multiplyByMInv(state,f,deltaU);
    if (u0.size()) state.updU() = u0 + deltaU;
    else state.updU() += deltaU;
    m_mbs.realize(state, Stage::Velocity);
}



//==============================================================================
//                       MY UNILATERAL CONSTRAINT SET
//==============================================================================


//------------------------ SELECT ACTIVE CONSTRAINTS ---------------------------
void MyUnilateralConstraintSet::
selectActiveConstraints(State& state, Real accuracy) const {

    // Find all the contacts and stiction elements that might be active based
    // on kinematics.
    MyElementSubset candidates;

    bool needRestart;
    do {
        //TODO: this (mis)use of accuracy needs to be revisited.
        findCandidateElements(state, accuracy, accuracy, candidates);

        // Evaluate accelerations and reaction forces and check if 
        // any of the active contacts are generating negative ("pulling") 
        // forces; if so, inactivate them.
        findActiveCandidates(state, candidates);

        Array_<Real> slipVels(candidates.m_friction.size());
        for (unsigned i=0; i < candidates.m_friction.size(); ++i) {
            const int which = candidates.m_friction[i];
            const MyFrictionElement& fric = getFrictionElement(which);
            slipVels[i] = fric.getActualSlipSpeed(state);
        }

        // Finally, project active constraints to the constraint manifolds.
        const Real tol = accuracy;
        SimTK_DEBUG1("projecting active constraints to tol=%g\n", tol);
        m_mbs.project(state, tol);

        // It is possible that the projection above took some marginally-sliding
        // friction and slowed it down enough to make it a stiction candidate.
        needRestart = false;
        for (unsigned i=0; i < candidates.m_sliding.size(); ++i) {
            const int which = candidates.m_sliding[i];
            const MyFrictionElement& fric = getFrictionElement(which);
            if (   fric.getActualSlipSpeed(state) <= accuracy
                || fric.calcSlipSpeedWitness(state) <= 0) 
            {
                SimTK_DEBUG3("***RESTART1** selectActiveConstraints() because "
                    "sliding velocity of friction %d is |v|=%g or witness=%g\n",
                    which, fric.getActualSlipSpeed(state),
                    fric.calcSlipSpeedWitness(state));
                needRestart = true;
                break;
            }
        }
        if (needRestart) continue;
        for (unsigned i=0; i < candidates.m_friction.size(); ++i) {
            const int which = candidates.m_friction[i];
            const MyFrictionElement& fric = getFrictionElement(which);
            if (fric.isSticking(state)) continue;
            if (slipVels[i] > accuracy
                && fric.getActualSlipSpeed(state) <= accuracy)
            {
                SimTK_DEBUG3("***RESTART2** selectActiveConstraints() because "
                    "sliding velocity of friction %d went from |v|=%g to |v|=%g\n",
                    which, slipVels[i], fric.getActualSlipSpeed(state));
                needRestart = true;
                break;
            }
        }

    } while (needRestart);
}




//-------------------------- FIND ACTIVE CANDIDATES ---------------------------
// Given a list of candidate unilateral constraints (contact and stiction),
// determine which ones are active in the least squares solution for the
// constraint multipliers. Candidates are those constraints that meet all 
// kinematic conditions -- for contacts, position and velocity errors less than
// tolerance; for stiction, sliding velocity less than tolerance. Also, any
// constraint that is currently active is a candidate, regardless of its
// kinematics.
//
// This method should be called only from within an event handler. For sliding
// friction it will have reset the "previous slip direction" to the current
// slip or impending slip direction, and converged the remembered normal force.
//
// Algorithm
// ---------
// We're given a set of candidates including contacts and stiction. If any are
// inactive, activate them.
// -- at this point all aerr==0, some ferr might be < 0
//
// loop
// - Realize(Acceleration) with the current active set
//     (iterate sliding forces until f=mu_d*N to tol)
// - Calculate ferr for active constraints, aerr for inactive
// - If all ferr>=0, aerr>=0 -> break loop
// - Check for aerr < 0 [tol?]. Shouldn't happen but if it does must turn on the
//     associated constraint for the worst violation, then -> continue loop
// - Find worst (most negative) offender:
//    contact offense  = fc < 0 ? fc : 0
//    stiction offense = mu_s*max(0, fc) - |fs|
// - Choose constraint to deactivate:
//     worst is a stiction constraint: choose it
//     worst is a contact constraint: if it has stiction, chose that
//                                    otherwise choose the contact constraint
// - Inactivate chosen constraint
//     (if stiction, record impending slip direction & N for stick->slide)
// end loop 
//
void MyUnilateralConstraintSet::
findActiveCandidates(State& s, const MyElementSubset& candidates) const
{
    const int ReviseNormalNIters = 6;
    showConstraintStatus(s, "ENTER findActiveCandidates()");
    if (candidates.m_contact.empty()) {
        // Can't be any friction either, if there are no contacts.
        SimTK_DEBUG("EXIT findActiveCandidates: no candidates.\n");
        m_mbs.realize(s, Stage::Acceleration);
        return;
    }

    SimTK_DEBUG3(
        "findActiveCandidates() for %d/%d/%d contact/stick/slip candidates ...\n",
        candidates.m_contact.size(), candidates.m_friction.size(),
        candidates.m_sliding.size());

    // Enable all candidate contact and stiction constraints if any are
    // currently disabled.
    enableConstraintSubset(candidates, true, s);

    int pass=0, nContactsDisabled=0, nStictionDisabled=0, nContactsReenabled=0;
    while (true) {
        ++pass; 
        SimTK_DEBUG1("\nfindActiveCandidates(): pass %d\n", pass);

        // Given an active set, evaluate multipliers and accelerations, and
        // converge sliding forces.
        m_mbs.realize(s, Stage::Acceleration);
        SimTK_DEBUG("findActiveCandidates(): CONVERGE NORMALS\n");
        convergeNormalForcesForSlidingContacts(s, SqrtEps);
 /*       for (int i=0; i < ReviseNormalNIters; ++i) {
            s.autoUpdateDiscreteVariables();
            s.invalidateAllCacheAtOrAbove(Stage::Dynamics);
            convergeNormalForcesForSlidingContacts(s, SqrtEps);
            m_mbs.realize(s, Stage::Acceleration);
        }*/
        m_mbs.realize(s, Stage::Acceleration);

        // Scan all candidate contacts to find the active one that has the
        // most negative normal force, and the inactive one that has the 
        // most negative acceleration error (hopefully none will).

        int worstActiveContact=-1; Real mostNegativeContactForce=0;
        int worstInactiveContact=-1; Real mostNegativeAerr=0;
        
        SimTK_DEBUG("analyzing contact constraints ...\n");
        for (unsigned i=0; i < candidates.m_contact.size(); ++i) {
            const int which = candidates.m_contact[i];
            SimTK_DEBUG1("  %d: ", which);
            const MyContactElement& cont = getContactElement(which);
            if (cont.isDisabled(s)) {
                const Real aerr = cont.getAerr(s);
                SimTK_DEBUG1("off aerr=%g\n", aerr);
                if (aerr < mostNegativeAerr) {
                    worstInactiveContact = which;
                    mostNegativeAerr = aerr;
                }
            } else {
                const Real f = cont.getForce(s);
                SimTK_DEBUG1("on f=%g\n", f);
                if (f < mostNegativeContactForce) {
                    worstActiveContact = which;
                    mostNegativeContactForce = f;
                }
            }
        }

        // This is bad and might cause cycling.
        if (mostNegativeAerr < 0) {
            SimTK_DEBUG2("  !!! Inactive contact %d violated, aerr=%g\n", 
                worstInactiveContact, mostNegativeAerr);
            const MyContactElement& cont = getContactElement(worstInactiveContact);
            //TODO -- must use a tolerance or prevent looping
            //++nContactsReenabled;
            //cont.enable(s);
            //continue;
        }

        SimTK_DEBUG("  All inactive contacts are satisfied.\n");

        #ifndef NDEBUG
        if (mostNegativeContactForce == 0)
            printf("  All active contacts are satisfied.\n");
        else 
            printf("  Active contact %d was worst violator with f=%g\n",
                worstActiveContact, mostNegativeContactForce);
        #endif

        int worstActiveStiction=-1; Real mostNegativeStictionCapacity=0;     
        SimTK_DEBUG("analyzing stiction constraints ...\n");
        for (unsigned i=0; i < candidates.m_friction.size(); ++i) {
            const int which = candidates.m_friction[i];
            SimTK_DEBUG1("  %d: ", which);
            const MyFrictionElement& fric = getFrictionElement(which);
            if (!fric.isSticking(s)) {
                SimTK_DEBUG("off\n");
                continue;
            }
            const Real mu_s = fric.getStaticFrictionCoef();
            const Real N = fric.getMasterNormalForce(s); // might be negative
            const Real fsmag = fric.getActualFrictionForce(s);
            const Real capacity = mu_s*std::max(N,Real(0)) - fsmag;
            SimTK_DEBUG2("on capacity=%g (N=%g)\n", capacity, N);

            if (capacity < mostNegativeStictionCapacity) {
                worstActiveStiction = which;
                mostNegativeStictionCapacity = capacity;
            }
        }

        #ifndef NDEBUG
        if (mostNegativeStictionCapacity == 0)
            printf("  All active stiction constraints are satisfied.\n");
        else 
            printf("  Active stiction %d was worst violator with capacity=%g\n",
                worstActiveStiction, mostNegativeStictionCapacity);
        #endif

        if (mostNegativeContactForce==0 && mostNegativeStictionCapacity==0) {
            SimTK_DEBUG("DONE. Current active set is a winner.\n");
            break;
        }

        if (mostNegativeStictionCapacity <= mostNegativeContactForce) {
            SimTK_DEBUG1("  Disable stiction %d\n", worstActiveStiction);
            MyFrictionElement& fric = updFrictionElement(worstActiveStiction);
            fric.recordImpendingSlipInfo(s);
            ++nStictionDisabled;
            fric.disableStiction(s);
            continue;
        }

        // A contact constraint was the worst violator. If that contact
        // constraint has an active stiction constraint, we have to disable
        // the stiction constraint first.
        SimTK_DEBUG1("  Contact %d was the worst violator.\n", worstActiveContact);
        const MyContactElement& cont = getContactElement(worstActiveContact);
        assert(!cont.isDisabled(s));

        if (cont.hasFrictionElement()) {
            MyFrictionElement& fric = cont.updFrictionElement();
            if (fric.isSticking(s)) {
                SimTK_DEBUG1("  ... but must disable stiction %d first.\n",
                    fric.getFrictionIndex());
                fric.recordImpendingSlipInfo(s);
                ++nStictionDisabled;
                fric.disableStiction(s);
                continue;
            }
        }

        SimTK_DEBUG1("  Disable contact %d\n", worstActiveContact); 
        ++nContactsDisabled;
        cont.disable(s);
        // Go back for another pass.
    }

    // Reset all the slip directions so that all slip->stick event witnesses 
    // will be positive when integration resumes.
    for (unsigned i=0; i < candidates.m_sliding.size(); ++i) {
        const int which = candidates.m_sliding[i];
        const MyFrictionElement& fric = getFrictionElement(which);
        if (!fric.isMasterActive(s)) continue;
        fric.updatePreviousSlipDirFromRecorded(s);
    }

    // Always leave at acceleration stage.
    m_mbs.realize(s, Stage::Acceleration);

    SimTK_DEBUG3("... Done; %d contacts, %d stictions broken; %d re-enabled.\n", 
        nContactsDisabled, nStictionDisabled, nContactsReenabled);

    showConstraintStatus(s, "EXIT findActiveCandidates()");
}

//---------------- CONVERGE NORMAL FORCES FOR SLIDING CONTACTS -----------------
namespace {
class NormalErrors : public Differentiator::JacobianFunction {
public:
    NormalErrors(const MyUnilateralConstraintSet& unis,
                 const Array_<int>&               sliders,
                 State&                           state) 
        : Differentiator::JacobianFunction(sliders.size(), sliders.size()), 
          m_unis(unis), m_state(state), m_sliders(sliders) { }

    void getEstimates(const State& s, Vector& estN) const {
        estN.resize(m_sliders.size());
        for (int i=0; i < (int)m_sliders.size(); ++i) {
            const MyFrictionElement& friction = 
                m_unis.getFrictionElement(m_sliders[i]);
            estN[i] = friction.getEstimatedNormalForce(s);
        }
    }

    // Alternate signature for convenience (extra copy).
    Vector getEstimates(const State& s) {
        Vector estN(m_sliders.size());
        getEstimates(s,estN);
        return estN;
    }

    void setEstimates(State& s, const Vector& estN) const {
        for (int i=0; i < (int)m_sliders.size(); ++i) {
            const MyFrictionElement& friction = 
                m_unis.getFrictionElement(m_sliders[i]);
            friction.setEstimatedNormalForce(m_state, estN[i]);
        }
    }

    void calcErrors(const State& s, Vector& errs) const {
        errs.resize(m_sliders.size());
        m_unis.getMultibodySystem().realize(s, Stage::Acceleration);
        for (int i=0; i < (int)m_sliders.size(); ++i) {
            const MyFrictionElement& friction = 
                m_unis.getFrictionElement(m_sliders[i]);
            const Real estN = friction.getEstimatedNormalForce(s);
            const Real N = friction.getMasterNormalForce(s);
            errs[i] = estN - N;
        }
    }

    // Alternate signature for convenience (extra copy).
    Vector calcErrors(const State& s) const {
        Vector errs(m_sliders.size()); calcErrors(s, errs);
        return errs;
    }

    // Must provide this pure virtual function.
    int f(const Vector& y, Vector& fy) const OVERRIDE_11 {
        assert(y.size() == m_sliders.size() && fy.size() == m_sliders.size());
        Vector prevEst(m_sliders.size());
        getEstimates(m_state, prevEst);
        setEstimates(m_state, y);
        calcErrors(m_state, fy);
        setEstimates(m_state, prevEst); // restore
        return 0;
    }
private:
    const MyUnilateralConstraintSet& m_unis;
    const Array_<int>&               m_sliders;
    State&                           m_state;
};
}


void MyUnilateralConstraintSet::
convergeNormalForcesForSlidingContacts(State& s, Real rtol) const {
    Array_<int> sliders;
    for (int i=0; i < getNumFrictionElements(); ++i) {
        const MyFrictionElement& friction = getFrictionElement(i);
        if (friction.isMasterActive(s) && !friction.isSticking(s))
            sliders.push_back(i);
    }
    SimTK_DEBUG1("convergeNormalForces: %d sliding contacts\n", sliders.size());
    if (sliders.empty()) return;

    NormalErrors normErrs(*this,sliders,s);
    Differentiator dnormErrs(normErrs);

    Vector estN(sliders.size()), err(sliders.size()), destN(sliders.size()),
           newEstN(sliders.size());
    normErrs.getEstimates(s, estN);
    normErrs.calcErrors(s, err);

    Real prevNorm = err.normInf();

    #ifndef NDEBUG
    std::cout << "====> Initial estN=" << estN << "-> err=" << err
         << " norm=" << prevNorm << std::endl;
    #endif

    if (prevNorm <= rtol)
        return;

    #ifndef NDEBUG
    std::cout << "\nNeed to converge sliding force and normal:\n";
    #endif

    Matrix J(sliders.size(), sliders.size());
    bool improved;
    const int MaxIters = 7;
    int iter = 0;
    do {
        ++iter;
        improved = false;
        dnormErrs.calcJacobian(estN, err, J, Differentiator::CentralDifference);
        FactorQTZ Jinv(J);
        #ifndef NDEBUG
        std::cout << "====> iter " << iter << " rank " << Jinv.getRank() 
             << " Jacobian\n" << J << std::endl;
        #endif  
        Jinv.solve(err, destN);
        Real fac = 2, norm;
        Real worseningAllowance = iter == 1 ? 1.2 : 1.01;
        do {
            fac /= 2;
            newEstN = estN - fac*destN;
            normErrs.setEstimates(s, newEstN);
            normErrs.calcErrors(s, err);
            norm = err.normInf();
            #ifndef NDEBUG
            printf("fac=%g new norm=%g", fac, norm);
            std::cout << " est=" << newEstN << std::endl;
            #endif  
        } while (norm > worseningAllowance*prevNorm);
        if (iter==1 || norm < prevNorm) {
            prevNorm = norm;
            estN = newEstN;
            improved = true;
        }
        #ifndef NDEBUG
        std::cout << "====> revised estN=" << estN << "-> err=" << err 
             << " norm=" << prevNorm << std::endl;
        #endif  
    } while (iter < MaxIters && improved && err.normInf() > rtol);
}


//-------------------------- SHOW CONSTRAINT STATUS ----------------------------
void MyUnilateralConstraintSet::
showConstraintStatus(const State& s, const String& place) const
{
#ifndef NDEBUG
    printf("\n%s: uni status @t=%.15g\n", place.c_str(), s.getTime());
    m_mbs.realize(s, Stage::Acceleration);
    for (int i=0; i < getNumContactElements(); ++i) {
        const MyContactElement& contact = getContactElement(i);
        const bool isActive = !contact.isDisabled(s);
        printf("  %6s %2d %s, h=%g dh=%g f=%g\n", 
                isActive?"ACTIVE":"off", i, contact.getContactType().c_str(), 
                contact.getPerr(s),contact.getVerr(s),
                isActive?contact.getForce(s):Zero);
    }
    for (int i=0; i < getNumFrictionElements(); ++i) {
        const MyFrictionElement& friction = getFrictionElement(i);
        if (!friction.isMasterActive(s))
            continue;
        const bool isSticking = friction.isSticking(s);
        printf("  %8s friction %2d, |v|=%g witness=%g\n", 
                isSticking?"STICKING":"sliding", i,
                friction.getActualSlipSpeed(s),
                isSticking?friction.calcStictionForceWitness(s)
                          :friction.calcSlipSpeedWitness(s));
        friction.writeFrictionInfo(s, "    ", std::cout);
    }
    printf("\n");
#endif
}


// This is the force handle's constructor; it just creates the force
// implementation object.
MySlidingFrictionForce::MySlidingFrictionForce
   (GeneralForceSubsystem&          forces,
    const MyPointContactFriction&   friction,
    Real                            vtol) 
:   Force::Custom(forces, new MySlidingFrictionForceImpl(forces,friction,vtol)) 
{}

Real MySlidingFrictionForce::getPrevN(const State& state) const 
{   return getImpl().getPrevN(state); }
void MySlidingFrictionForce::setPrevN(State& state, Real N) const 
{   getImpl().setPrevN(state,N); }
Vec2 MySlidingFrictionForce::getPrevSlipDir(const State& state) const 
{   return getImpl().getPrevSlipDir(state); }
bool MySlidingFrictionForce::hasPrevSlipDir(const State& state) const 
{   return !getImpl().getPrevSlipDir(state).isNaN(); }
void MySlidingFrictionForce::
setPrevSlipDir(State& state, const Vec2& slipDir) const
{   getImpl().setPrevSlipDir(state, slipDir); }
Real MySlidingFrictionForce::calcSlidingForceMagnitude(const State& state) const 
{   return getImpl().calcSlidingForceMagnitude(state); }
Vec2 MySlidingFrictionForce::calcSlidingForce(const State& state) const 
{   return getImpl().calcSlidingForce(state); }


const MySlidingFrictionForceImpl& MySlidingFrictionForce::
getImpl() const {
    return dynamic_cast<const MySlidingFrictionForceImpl&>(getImplementation()); 
}


//    --------------------------------
//    perr = ~p_BS*n - h
//    --------------------------------
void UniPointInPlaneImpl::calcPositionErrors      
   (const State&                                    s,      // Stage::Time
    const Array_<Transform,ConstrainedBodyIndex>&   allX_AB, 
    const Array_<Real,     ConstrainedQIndex>&      constrainedQ,
    Array_<Real>&                                   perr)   // mp of these
    const 
{
    assert(allX_AB.size()==2 && constrainedQ.size()==0 && perr.size() == 1);

    const Vec3       p_AS = findStationLocation(allX_AB, followerBody, 
                                                defaultFollowerPoint);
    const Transform& X_AB = getBodyTransform(allX_AB, planeBody);
    const Vec3       p_BC = ~X_AB * p_AS; // shift to B origin, reexpress in B;
                                     // C is material pt of B coincident with S

    // We'll calculate this scalar using B-frame vectors, but any frame would 
    // have done.
    perr[0] = dot(p_BC, defaultPlaneNormal) - defaultPlaneHeight;
}

//    --------------------------------
//    verr = ~v_CS_A*n
//    --------------------------------
void UniPointInPlaneImpl::calcPositionDotErrors      
   (const State&                                    s,      // Stage::Position
    const Array_<SpatialVec,ConstrainedBodyIndex>&  V_AB, 
    const Array_<Real,      ConstrainedQIndex>&     constrainedQDot,
    Array_<Real>&                                   pverr)  // mp of these
    const
{
    assert(V_AB.size()==2 && constrainedQDot.size()==0 && pverr.size() == 1);
    //TODO: should be able to get p info from State

    const Vec3       p_AS = findStationLocationFromState(s, followerBody, 
                                                         defaultFollowerPoint);
    const Transform& X_AB = getBodyTransformFromState(s, planeBody);
    const Vec3       p_BC = ~X_AB * p_AS; // shift to B origin, reexpress in B;
                                   // C is material point of B coincident with S
    const UnitVec3   n_A  = X_AB.R() * defaultPlaneNormal;

    const Vec3       v_AS = findStationVelocity(s, V_AB, followerBody, 
                                                defaultFollowerPoint);
    const Vec3       v_AC = findStationVelocity(s, V_AB, planeBody, p_BC);

    // Calculate this scalar using A-frame vectors.
    pverr[0] = dot( v_AS-v_AC, n_A );
}

//    -------------------------------------
//    aerr = ~(a_CS_A - 2 w_AB X v_CS_A) * n
//    -------------------------------------
void UniPointInPlaneImpl::calcPositionDotDotErrors      
   (const State&                                    s,      // Stage::Velocity
    const Array_<SpatialVec,ConstrainedBodyIndex>&  A_AB, 
    const Array_<Real,      ConstrainedQIndex>&     constrainedQDotDot,
    Array_<Real>&                                   paerr)  // mp of these
    const
{
    assert(A_AB.size()==2 && constrainedQDotDot.size()==0 && paerr.size() == 1);
    //TODO: should be able to get p and v info from State
    const Vec3       p_AS = findStationLocationFromState(s, followerBody, 
                                                         defaultFollowerPoint);
    const Transform& X_AB = getBodyTransformFromState(s, planeBody);
    const Vec3       p_BC = ~X_AB * p_AS; // shift to B origin, reexpress in B;
                                   // C is material point of B coincident with S
    const UnitVec3   n_A  = X_AB.R() * defaultPlaneNormal;

    const Vec3&      w_AB = getBodyAngularVelocityFromState(s, planeBody);
    const Vec3       v_AS = findStationVelocityFromState(s, followerBody, 
                                                         defaultFollowerPoint);
    const Vec3       v_AC = findStationVelocityFromState(s, planeBody, p_BC);

    const Vec3       a_AS = findStationAcceleration(s, A_AB, followerBody, 
                                                    defaultFollowerPoint);;
    const Vec3       a_AC = findStationAcceleration(s, A_AB, planeBody, p_BC);

    paerr[0] = dot( (a_AS-a_AC) - 2.*w_AB % (v_AS-v_AC), n_A );
}
