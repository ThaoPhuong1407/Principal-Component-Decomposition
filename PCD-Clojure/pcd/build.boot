(def project 'pcd)
(def version "0.1.0-SNAPSHOT")

(set-env! :resource-paths #{"resources" "src"}
          :source-paths   #{"test"}
          :dependencies   '[[org.clojure/clojure "1.10.0"]
                            ;; [adzerk/boot-test "RELEASE" :scope "test"]
                            [org.clojure/data.csv "0.1.3"]
                            [net.mikera/core.matrix "0.62.0"]
                            [incanter "1.9.3"]
                            [clj-python/libpython-clj "2.00-beta-15"]
                            #_[cnuernber/libpython-clj "1.36"]])

(task-options!
 aot {:namespace   #{'pcd.core}}
 pom {:project     project
      :version     version
      :description "FIXME: write description"
      :url         "http://example/FIXME"
      :scm         {:url "https://github.com/yourname/pcd"}
      :license     {"Eclipse Public License"
                    "http://www.eclipse.org/legal/epl-v10.html"}}
 repl {:init-ns    'pcd.core}
 jar {:main        'pcd.core
      :file        (str "pcd-" version "-standalone.jar")})

(deftask build
  "Build the project locally as a JAR."
  [d dir PATH #{str} "the set of directories to write to (target)."]
  (let [dir (if (seq dir) dir #{"target"})]
    (comp (aot) (pom) (uber) (jar) (target :dir dir))))

(deftask run
  "Run the project."
  [a args ARG [str] "the arguments for the application."]
  (with-pass-thru fs
    (require '[pcd.core :as app])
    (apply (resolve 'app/-main) args)))

;; (require '[adzerk.boot-test :refer [test]])
