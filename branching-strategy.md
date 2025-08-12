Branching Strategy
========

## Branch from 'development'

1. **Branching from `development`:**
   - Create feature branches from the `development` branch for your work. Do this by using the following commands:
     ```bash
     git checkout development
     git pull origin development
     git checkout -b feature/your-feature-name
     ```

2. **Creating Pull Requests:**
   - Once you have completed your work on a feature branch, you should push your branch to the remote repository:
     ```bash
     git push origin feature/your-feature-name
     ```
   - Then, you can create a pull request (PR) in Azure DevOps:
     - Go to the **Pull Requests** section in your repository.
     - Click on **New pull request**.
     - Set the **Source branch** to your feature branch and the **Target branch** to `development`.
     - Fill in the PR details and submit it for review.

## Merge into `main`

1. **Review and Merge PRs:**
   - Once the PR against the `development` branch is approved and passes any policies, it can be merged into `development`.
   - After testing and validating the changes in `development`, you can create a PR to merge `development` into `main` when ready for deployment.

2. **Merge `development` into `main`:**
   - Create a new PR with the **Source branch** set to `development` and the **Target branch** set to `main`.
   - Review and merge the PR once it’s approved.

### Step 5: Regularly Sync `development` with `main`

- To keep the `development` branch up to date with any changes made directly to `main`, regularly merge `main` into `development`:
  ```bash
  git checkout development
  git pull origin main
  ```
