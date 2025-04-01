package modi;

import msutil.ProtCutter;

import java.util.concurrent.LinkedBlockingQueue;
import java.util.concurrent.ThreadFactory;
import java.util.concurrent.ThreadPoolExecutor;
import java.util.concurrent.TimeUnit;
import java.util.concurrent.atomic.AtomicInteger;

public class Mutables {
    public static ProtCutter protease = ProtCutter.getCutter("Trypsin");
    public static double fragmentTolerance = 0.6; // never change after set param invocation
    public static PTMDB variableModifications;
    public static PTMDB fixedModifications;
    public static double[] reporterMassOfIsobaricTag = null;
    public static double massToleranceForDenovo = 0.3;

    public static int    maxNoOfC13 = 0;
    public static double precursorAccuracy = 0.5;
    public static double precursorTolerance = 0.5;
    public static double gapTolerance = 0.6;
    public static double gapAccuracy = 1.6;
    public static double nonModifiedDelta = massToleranceForDenovo;






    public static boolean fEqual(double v1, double v2) {
        return Math.abs(v1 - v2) <= Mutables.fragmentTolerance;
    }


    @SuppressWarnings("BooleanMethodIsAlwaysInverted")
    public static boolean isInModifiedRange(double v) {
        if (Constants.minModifiedMass - Mutables.gapTolerance < v && v < Constants.maxModifiedMass + Mutables.gapTolerance)
            return true;
        else return Math.abs(v) <= Mutables.gapTolerance;
    }

    public static boolean isWithinTolerance(double calc, double obsv, double tol) {

        if (Constants.minNoOfC13 == 0 && Mutables.maxNoOfC13 == 0) {
            return !(Math.abs(calc - obsv) > tol);
        } else {
            double tempError = obsv - calc;
            int isoerr = Constants.round(tempError / Constants.IsotopeSpace);
            if (isoerr < Constants.minNoOfC13 || Mutables.maxNoOfC13 < isoerr) return false;
            return !(Math.abs(tempError - isoerr * Constants.IsotopeSpace) > Mutables.precursorAccuracy);
        }
    }

    public static boolean isWithinAccuracy(double err) {
        if (Mutables.gapAccuracy > 0.5) return true;
        int isoerr = Constants.round(err / Constants.IsotopeSpace);
        return !(Math.abs(err - isoerr * Constants.IsotopeSpace) > Mutables.gapAccuracy);
    }

    public static void printGlobalMutables() {
        System.out.println("*****************************");
        System.out.println("*****************************");
        System.out.println("maxNoOfC13: " + maxNoOfC13); //analyze
        System.out.println("precursorAccuracy: " + precursorAccuracy);
        System.out.println("precursorTolerance: " + precursorTolerance);
        System.out.println("gapTolerance: " + gapTolerance);
        System.out.println("gapAccuracy: " + gapAccuracy);
        System.out.println("nonModifiedDelta: " + nonModifiedDelta);
        System.out.println("*****************************");
        System.out.println("*****************************");
    }
/*임시비활(하지만 곧 적용할 코드임)




    //reference(actual practical usage example)
    public class CustomWorker extends Thread {

        public final ThreadLocal<String> CUSTOM_THREAD_LOCAL = new InheritableThreadLocal<>() {
            @Override
            protected String childValue(String parentValue) {
                return "child-" + parentValue;
            }
        };

        // Runnable을 받는 생성자
        public CustomWorker(Runnable target) {
            super(target);
        }

        @Override
        public void run() {
            // ThreadLocal 값 설정
            CUSTOM_THREAD_LOCAL.set("custom-thread-value");
            System.out.println(getName() + " - ThreadLocal value: " + CUSTOM_THREAD_LOCAL.get());

            // 실제 작업 실행 (Thread 클래스의 run 메서드는 target.run()을 호출)
            super.run();

            // 작업 완료 후 ThreadLocal 정리
            CUSTOM_THREAD_LOCAL.remove();
        }
    }


    class CustomWorkerFactory implements ThreadFactory {
        private static final AtomicInteger poolNumber = new AtomicInteger(1);
        private final AtomicInteger threadNumber = new AtomicInteger(1);
        private final String namePrefix;

        public CustomWorkerFactory() {
            namePrefix = "custom-pool-" + poolNumber.getAndIncrement() + "-thread-";
        }

        @Override
        public Thread newThread(Runnable r) {
            CustomWorker thread = new CustomWorker(r);
            thread.setName(namePrefix + threadNumber.getAndIncrement());
            return thread;
        }
    }

    void worker_pool_multi_tasks_perform_example(){

        ThreadPoolExecutor executor = new ThreadPoolExecutor(
                10,                          // 코어 스레드 수
                10,                         // 최대 스레드 수
                0L,                        // 유휴 스레드 유지 시간
                TimeUnit.MILLISECONDS,           // 시간 단위
                new LinkedBlockingQueue<>(10000),  // 작업 큐
                new CustomWorkerFactory(),  // 커스텀 스레드 팩토리
                new ThreadPoolExecutor.CallerRunsPolicy()  // 거부 정책
        );

        // 여러 작업 제출(시뮬레이션)
        for (int i = 0; i < 20; i++) {
            final int taskId = i;
            executor.execute(() -> {
                String threadName = Thread.currentThread().getName();
                //String localValue = CustomWorker.CUSTOM_THREAD_LOCAL.get(); // no race condition 일까 과연? 각각 로컬 값 가질까? 완벽할까? 이 줄 테스트.

                System.out.println("작업 " + taskId + " 실행 중 - 스레드: " + threadName);
                //System.out.println("작업 " + taskId + " - ThreadLocal 값: " + localValue);

                // 특정 함수 내부라고 가정하고, currentThread().. 호출.. 및 스레드로컬 컨텍스트 접근 해보기
                ((CustomWorker)Thread.currentThread()).CUSTOM_THREAD_LOCAL.set("");  // read and write 다 테스트 필요

                try {
                    // 의도적인 시뮬레이션
                    Thread.sleep(1000);
                } catch (InterruptedException e) {
                    Thread.currentThread().interrupt();
                }


                System.out.println("[Task " + taskId + " completed]");
            });
        }

        executor.shutdown();
        try {
            // 모든 작업이 완료될 때까지 최대 30초 대기
            if (!executor.awaitTermination(30, TimeUnit.SECONDS)) {
                System.out.println("일부 작업이 완료되지 않았습니다. 강제 종료합니다.");
                executor.shutdownNow();
            }
        } catch (InterruptedException e) {
            //System.out.println("종료 대기 중 인터럽트 발생");
            executor.shutdownNow();
            Thread.currentThread().interrupt();
        }

        System.out.println("모든 작업 완료");
    }






    */
     //
     //
     //
     //


    //reference(basic)
    /*

    class CustomWorker extends Thread {
        public static final ThreadLocal<String> INHERITABLE_THREAD_LOCAL = new InheritableThreadLocal<>() {
            @Override
            public String childValue(String parentValue) {
                return "child value";
            }
        };


        public CustomWorker(Runnable target) {
            super(target);
        }

        @Override
        public void run() {
            INHERITABLE_THREAD_LOCAL.set("CustomWorker");
            System.out.println("CustomWorker Thread : " + Thread.currentThread());
            System.out.println("CustomWorker ThreadLocal Variable : " + INHERITABLE_THREAD_LOCAL.get());

            // 상위 클래스의 run 메서드 호출 (제출된 작업 실행)
            super.run();  // <- 이 부분만 조금 특별히 다름. caution 주의 필요..


        }

    }


    class SimpleParentThreadFactory implements ThreadFactory {
        @Override
        public Thread newThread(Runnable r) {
            return new Parent(r);
        }
    }






     */

}
