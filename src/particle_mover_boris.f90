module m_particle_mover_boris
    use m_vector, only: cross
    use oh_type, only: oh_particle
    use m_particle_mover, only: t_ParticleMover
    implicit none

    private
    public t_ParticleMoverBoris, new_ParticleMoverBoris

    type, extends(t_ParticleMover) :: t_ParticleMoverBoris
    contains
        procedure :: move => boris_move
    end type

contains

    function new_ParticleMoverBoris() result(obj)
        type(t_ParticleMoverBoris) :: obj
    end function

    subroutine boris_move(self, particle, qm, eb, dt)
        class(t_ParticleMoverBoris), intent(in) :: self
        type(oh_particle), intent(inout) :: particle
        double precision, intent(in) :: qm
        double precision, intent(in) :: eb(6)
        double precision :: dt

        block
            double precision :: ef(3), bf(3)
            double precision :: uold(3), unew(3)
            double precision :: upm(3), upa(3), upp(3)
            double precision :: s(3), t(3)

            double precision :: dt2

            ! Velocity update
            dt2 = dt*0.5d0

            ef(:) = eb(1:3)
            bf(:) = eb(4:6)

            t(:) = qm*bf(:)*dt2
            s(:) = 2*t(:)/(1 + sum(t(:)*t(:)))

            uold(:) = [particle%vx, particle%vy, particle%vz]

            upm(:) = uold(:) + qm*ef(:)*dt2

            upa(:) = upm(:) + cross(upm, t)
            upp(:) = upm(:) + cross(upa, s)

            unew(:) = upp(:) + qm*ef(:)*dt2

            particle%vx = unew(1)
            particle%vy = unew(2)
            particle%vz = unew(3)
        end block

        ! Position update
        particle%x = particle%x + particle%vx*dt
        particle%y = particle%y + particle%vy*dt
        particle%z = particle%z + particle%vz*dt

    end subroutine

end module
