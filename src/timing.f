      subroutine timing(t)

c Returns the time in seconds.
c This call is a system-dependent call

      double precision t

      t=0.d0

      t=dble(mclock())/1.d6
c      call cpu_time(t)

      return

      end
